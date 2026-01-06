(defvar *counter* 0)

(defparameter *default-contraction-rules*
  '((W . contract-W)
    (VARPHI . contract-varphi))
  "Default contraction rules for the NLSM.")

(defparameter *default-fields*
  '(W VARPHI)
  "Default field types to look for in the Keldysh electronic model.")

(defvar *contraction-rules*
  "Dynamic configuration map. 
   format: ((TYPE . FUNCTION-NAME) ...)
   Example: ((W . contract-W) (varphi . contract-varphi))")

;; ============================================================
;; 1. NORMALIZATION & UTILS
;; ============================================================

(defun expand-trace-arg (term)
  "Expands an exponent expression into a list of repeated terms.
   
   Positive Exponent:
   Input: (expt X 3) -> Output: (X X X)
   
   Negative Exponent:
   Input: (expt X -2) -> Output: ((/ 1 (* X X)))
   
   If input is not an exponent, returns it as a single-item list."
  (if (and (consp term) (eq (first term) 'expt))
      (let ((base (second term))
            (exponent (third term)))
        (if (minusp exponent)
            (let* ((n (abs exponent))
                   (repeated (loop repeat n collect (copy-tree base))))
              (list (list '/ 1 (cons '* repeated))))
            (loop repeat exponent collect (copy-tree base))))
      (list term)))

(defun normalize (expr)
  "Walks the entire expression tree to handle pre-processing tasks.
   Currently, it expands powers inside traces (Tr[X^2] -> Tr[X X])."
  (cond
    ((atom expr) expr)
    ((eq (first expr) 'tr)
     (cons 'tr (mapcan #'expand-trace-arg (cdr expr))))
    (t (cons (first expr) (mapcar #'normalize (cdr expr))))))

(defun remove-nth (n list)
  "Removes the item at index N from the list (non-destructively)."
  (if (zerop n) (cdr list)
      (cons (car list) (remove-nth (1- n) (cdr list)))))

(defun deep-copy (expr) 
  "Creates a fresh deep copy of the expression tree."
  (copy-tree expr))

;; ============================================================
;; 2. TAGGING
;; ============================================================

(defun get-tag-id (tag)
  "Retrieves the unique numeric ID from a tagged field.
   Input: (:tagged ID TYPE CONTENT) -> Output: ID (integer)"
  (second tag))

(defun get-tag-type (tag)
  "Retrieves the field type symbol from a tagged field.
   Input: (:tagged ID TYPE CONTENT) -> Output: TYPE (symbol)"
  (third tag))

(defun get-tag-content (tag)
  "Retrieves the actual expression content from a tagged field.
   Input: (:tagged ID TYPE CONTENT) -> Output: CONTENT (list)"
  (fourth tag))

(defun tag-fields (expr target-fields)
  "Recursively walks the expression tree. If it encounters a symbol belonging
   to the list TARGET-FIELDS, it wraps it in a tagged structure including its type.
   
   Input: 
     expr: The expression tree.
     target-fields: A list of symbols to target (e.g., '(W VARPHI))
   
   Output: 
     The expression tree with targets replaced by (:tagged ID TYPE CONTENT)."
  (cond
    ((and (consp expr) (member (first expr) target-fields))
     (incf *counter*)
     (list :tagged *counter* (first expr) (copy-tree expr)))
    ((consp expr)
     (cons (tag-fields (car expr) target-fields)
           (tag-fields (cdr expr) target-fields)))
    (t expr)))

(defmacro with-tagging (&body body)
  "Sets up the tagging environment (resets the counter)."
  `(let ((*counter* 0)) ,@body))

(defun extract-tags (expr)
  "Returns a flat list of all (:tagged ...) objects found in the tree."
  (cond
    ((atom expr) nil)
    ((eq (first expr) :tagged) (list expr))
    (t (mapcan #'extract-tags expr))))

;; ============================================================
;; 3. PAIRING LOGIC
;; ============================================================

(defun group-tags-by-type (tags)
  "Sorts a flat list of tags into buckets based on their type.
   Input: (tag1 tag2 ...)
   Output: Alist ((Type1 . (tags...)) (Type2 . (tags...)))"
  (let ((groups nil))
    (dolist (tag tags)
      (let* ((type (get-tag-type tag))
             (entry (assoc type groups)))
        (if entry
            (push tag (cdr entry))
            (push (list type tag) groups))))
    groups))

(defun cartesian-product (lists-of-pairings)
  "Recursively combines independent pairing lists into global scenarios."
  (cond
    ((null lists-of-pairings) '(nil))
    (t (let ((current-scenarios (first lists-of-pairings))
             (rest-scenarios (cartesian-product (rest lists-of-pairings))))
         (loop for curr in current-scenarios
               nconc (loop for r in rest-scenarios
                           collect (append curr r)))))))

(defun generate-pairs (items)
  "Generates ALL possible Wick contraction pairings for the list of items.
   Recursive logic: Pick head, pair with every other item, recurse on remainder."
  (cond
    ((null items) '(nil))
    ((oddp (length items)) (error "Odd number of fields."))
    (t
     (let ((head (first items))
           (rest (cdr items)))
       ;; Iterate through all possible partners for 'head'
       (loop for i from 0 below (length rest)
             for partner in rest
             for remaining = (remove-nth i rest)
             ;; Combine this pair with all valid pairings of the remaining items
             nconc (loop for sub-pairing in (generate-pairs remaining)
                         collect (cons (list head partner) sub-pairing)))))))

(defun generate-stratified-pairing (items)
  "The Main Pairing Function.
   Input: items (A flat list of tagged fields).
   Output: List of all valid global pairing."
  (let* ((groups (group-tags-by-type items))
         (per-type-scenarios 
          (loop for entry in groups
                collect (generate-pairs (cdr entry)))))
    (cartesian-product per-type-scenarios)))

;; ============================================================
;; 4. FIELD CONTRACTION LOGIC
;; ============================================================

(defun apply-op (op list-of-terms)
  "Wraps a list of terms in an operator (PARA/PERP).
   Handles the edge case where the list is empty -> returns Identity."
  (if (null list-of-terms) 
      (list (list op 'ID)) 
      (list (cons op list-of-terms))))

(defun find-trace-containing (expr id)
  "Locates the specific (TR ...) list object in the skeleton that contains
   the tag with the given ID. Returns the list object itself."
  (cond
    ((atom expr) nil)
    ((eq (first expr) 'tr)
     (if (member-if (lambda (x) (and (consp x) (eq (first x) :tagged) (= (second x) id))) expr)
         expr nil))
    (t (some (lambda (x) (find-trace-containing x id)) expr))))

(defun get-content-without-tr (trace-obj) 
  "Helper to strip the 'TR symbol from a trace list."
  (cdr trace-obj))

(defun get-single-coefficient (trace-content id)
  "for Different-Trace contraction (Tr[A W]).
   Cyclically rotates the trace content so that W (id) is at the end,
   then returns everything else (A)."
  (let ((pos (position-if
	      (lambda (x) (and (consp x) (eq (first x) :tagged) (= (second x) id))) trace-content)))
    (append (subseq trace-content (1+ pos)) (subseq trace-content 0 pos))))

(defun get-splitting-coefficients (trace-content id1 id2)
  "for Same-Trace contraction (Tr[A W1 B W2]).
   Returns two values: the list A (between W1 and W2) and the list B (wrapping around)."
  (let ((pos1 (position-if (lambda (x) (and (consp x) (eq (first x) :tagged) (= (second x) id1))) trace-content))
        (pos2 (position-if (lambda (x) (and (consp x) (eq (first x) :tagged) (= (second x) id2))) trace-content)))
    (if (< pos1 pos2)
        (values (subseq trace-content (1+ pos1) pos2)
                (append (subseq trace-content (1+ pos2)) (subseq trace-content 0 pos1)))
        (values (append (subseq trace-content (1+ pos1)) (subseq trace-content 0 pos2))
                (subseq trace-content (1+ pos2) pos1)))))

(defun contract-W (skeleton tag1 tag2)
  "The specific physics logic for Standard Fermionic Fields (W).
   Handles: Derivatives extraction, Propagator construction, Trace Splitting/Merging.
   Input: skeleton (tree), tag1, tag2 -> Output: S-Expression of the contraction."
  (let* ((id1 (get-tag-id tag1)) (id2 (get-tag-id tag2))
         (content1 (get-tag-content tag1)) 
         (content2 (get-tag-content tag2))
         
         ;; --- 1. Extract Derivatives ---
         (r1 (second content1)) (r2 (second content2))
         (all-derivs (append (cddr content1) (cddr content2)))
         
         ;; --- 2. Construct Propagator ---
         (propagator (list* 'PI (list '- r1 r2) all-derivs))
         
         ;; --- 3. Locate Traces ---
         (trace1 (find-trace-containing skeleton id1))
         (trace2 (find-trace-containing skeleton id2)))

    (cond
      ((or (null trace1) (null trace2)) (error "Tags missing in skeleton."))
      
      ;; Case A: SPLITTING (Same Trace)
      ((eq trace1 trace2)
       (multiple-value-bind (A-list B-list) 
           (get-splitting-coefficients (get-content-without-tr trace1) id1 id2)
         (let ((A-para (apply-op 'PARA (deep-copy A-list)))
               (B-para (apply-op 'PARA (deep-copy B-list))))
           `(* (/ (- t-coupling) 4)
               (- (* (tr ,@A-para ,propagator) (tr ,@B-para))
                  (* (tr ,@(deep-copy A-para) sigma3 ,(deep-copy propagator)) 
                     (tr ,@(deep-copy B-para) sigma3)))))))

      ;; Case B: MERGING (Different Traces)
      (t 
       (let* ((A-list (get-single-coefficient (get-content-without-tr trace1) id1))
              (B-list (get-single-coefficient (get-content-without-tr trace2) id2))
              (A-perp (apply-op 'PERP (deep-copy A-list)))
              (B-perp (apply-op 'PERP (deep-copy B-list))))
         `(* (/ (- t-coupling) 2) (tr ,@A-perp ,propagator ,@B-perp)))))))

(defun replace-tag-in-expr (expr target-id replacement-list)
  "Recursively walks the expression tree.
   When it finds a Tag with ID = target-id, it replaces the Tag 
   with the contents of replacement-list (spliced in)."
  (cond
    ;; Case: Found the target tag -> Return Splice Instruction
    ((and (consp expr) (eq (first expr) :tagged) (= (second expr) target-id))
     (cons :splice replacement-list))

    ;; Case: List - Handle Splicing
    ((consp expr)
     (loop for item in expr
           for processed = (replace-tag-in-expr item target-id replacement-list)
           if (and (consp processed) (eq (first processed) :splice))
             append (cdr processed)
           else collect processed))

    (t expr)))

(defun contract-varphi (skeleton tag1 tag2)
  "Contraction logic for Varphi fields.
   Generates 3 branches (VK, VA, VR).
   UPDATED: Moves the scalar potentials (VK/VA/VR) OUT of the trace 
   and groups them with the prefactor."
  (let* ((id1 (get-tag-id tag1)) (id2 (get-tag-id tag2))
         (content1 (get-tag-content tag1)) 
         (content2 (get-tag-content tag2))
         (x (second content1))
         (y (second content2))
         (diff (list '- x y))
         (prefactor '(/ i 2))
         
         ;; NORMALIZE SKELETON:
         (skel-list (if (eq (first skeleton) 'tr)
                        (list skeleton)
                        skeleton)))

    (list '+ 
          ;; --- Branch 1: VK (Scalar: VK, Matrix: 1, 1) ---
          (let ((skel (deep-copy skel-list)))
            (setf skel (replace-tag-in-expr skel id1 (list 1))) ;; Just Matrix
            (setf skel (replace-tag-in-expr skel id2 (list 1))) ;; Just Matrix
            `(* ,prefactor ,(list 'VK diff) ,@skel))

          ;; --- Branch 2: VA (Scalar: VA, Matrix: SIGMA1, 1) ---
          (let ((skel (deep-copy skel-list)))
            (setf skel (replace-tag-in-expr skel id1 (list 'SIGMA1)))
            (setf skel (replace-tag-in-expr skel id2 (list 1)))
            `(* ,prefactor ,(list 'VA diff) ,@skel))

          ;; --- Branch 3: VR (Scalar: VR, Matrix: 1, SIGMA1) ---
          (let ((skel (deep-copy skel-list)))
            (setf skel (replace-tag-in-expr skel id1 (list 1)))
            (setf skel (replace-tag-in-expr skel id2 (list 'SIGMA1)))
            `(* ,prefactor ,(list 'VR diff) ,@skel)))))


;; ============================================================
;; 5. PHYSICS RULES DISPATCHER
;; ============================================================

(defun apply-contraction-rules (skeleton pairing)
  "The Dispatcher. Looks up the correct contraction function based on Tag Type."
  (let* ((tag1 (first pairing))
         (tag2 (second pairing))
         (type (get-tag-type tag1)) ;; We assume tag1 and tag2 have same type
         (rule-entry (assoc type *contraction-rules*)))
    
    (unless rule-entry
      (error "No contraction rule found for field type: ~A" type))
    
    (funcall (cdr rule-entry) skeleton tag1 tag2)))

;; ============================================================
;; 6. EXPANSION ENGINE
;; ============================================================

(defun has-tagged-field-p (expr)
  "Checks if an expression tree contains any tagged fields."
  (cond ((atom expr) nil)
        ((eq (first expr) :tagged) t)
        (t (some #'has-tagged-field-p expr))))

(defun distribute-linear-op (op args)
  "Distributes a linear operator (Tr, Para, Perp) over a Sum.
   Also handles flattening products inside traces if necessary."
  (let ((pos-of-sum (position-if (lambda (x) (and (consp x) (eq (first x) '+))) args)))
    (if (null pos-of-sum) (cons op args)
        (let ((sum-terms (cdr (nth pos-of-sum args)))
              (prefix (subseq args 0 pos-of-sum))
              (suffix (subseq args (1+ pos-of-sum))))
          (cons '+ (loop for term in sum-terms collect
                         (let* ((middle (if (and (eq op 'tr) (consp term) (eq (first term) '*))
                                            (cdr term) (list term))))
                           (expand-expression (cons op (append (deep-copy prefix) middle (deep-copy suffix)))))))))))

(defun make-op-term (op content)
  "Creates (OP content...). 
   Fixes empty content: PARA -> 'ID, PERP -> 0."
  (if (null content) 
      (if (eq op 'PARA) 'ID 0)
      (cons op content)))

(defun expand-operator (op terms)
  "Expands (OP A W B) based on the parity of W (assumed Odd/Perp).
   Recursive expansion ensures multiple Ws are handled correctly."
  (if (not (has-tagged-field-p terms)) 
      (cons op (deep-copy terms))
      
      (let* ((pos (position-if (lambda (x) (and (consp x) (eq (first x) :tagged))) terms))
             (W (nth pos terms))
             (L (subseq terms 0 pos))
             (R (subseq terms (1+ pos))))

        ;; W is inherently PERP (Odd).
        ;; PARA (Even) -> L*R must be Odd (Mixed Parity).
        ;; PERP (Odd)  -> L*R must be Even (Same Parity).
        (if (eq op 'PARA)
            `(+ (* ,(make-op-term 'PARA (deep-copy L)) 
                   ,(deep-copy W) 
                   ,(expand-expression (make-op-term 'PERP (deep-copy R)))) 
                (* ,(make-op-term 'PERP (deep-copy L)) 
                   ,(deep-copy W) 
                   ,(expand-expression (make-op-term 'PARA (deep-copy R)))))

            `(+ (* ,(make-op-term 'PARA (deep-copy L)) 
                   ,(deep-copy W) 
                   ,(expand-expression (make-op-term 'PARA (deep-copy R))))
                (* ,(make-op-term 'PERP (deep-copy L)) 
                   ,(deep-copy W) 
                   ,(expand-expression (make-op-term 'PERP (deep-copy R)))))))))

(defun expand-expression (expr)
  "The Main Expansion Routine.
   Handles: Minus signs, Flattening Sums, Distributing operators, Para/Perp expansion."
  (cond
    ((atom expr) expr)
    ((eq (first expr) :tagged) expr)
    
    ;; Handle Minus: (- A) -> (* -1 A)
    ((eq (first expr) '-)
     (let ((args (cdr expr)))
       (if (= (length args) 1) (expand-expression `(* -1 ,(first args)))
           (expand-expression (cons '+ (cons (first args) (mapcar (lambda (x) `(* -1 ,x)) (rest args))))))))

    ;; Handle Sums: Flatten nested additions
    ((eq (first expr) '+)
     (let ((expanded (mapcar #'expand-expression (cdr expr))))
       (cons '+ (loop for term in expanded 
                      if (and (consp term) (eq (first term) '+)) 
                      append (cdr term) 
                      else collect term))))

    ;; Handle Products & Traces: Distribute linearity
    ((or (eq (first expr) '*) (eq (first expr) 'tr))
     (let ((args (mapcar #'expand-expression (cdr expr))))
       (distribute-linear-op (first expr) 
                             (if (eq (first expr) '*) 
                                 ;; Flatten nested products
                                 (loop for a in args if (and (consp a) (eq (first a) '*)) append (cdr a) else collect a)
                                 args))))

    ;; Handle Para/perp operator: Trigger expansion if W is inside
    ((or (eq (first expr) 'PARA) (eq (first expr) 'PERP))
     (if (has-tagged-field-p expr)
         (expand-expression (expand-operator (first expr) (cdr expr)))
         (distribute-linear-op (first expr) (mapcar #'expand-expression (cdr expr)))))

    (t expr)))

;; ============================================================
;; 7. ORCHESTRATOR
;; ============================================================

(defun contains-tag-id-p (expr id)
  "Helper to find if an expression subtree contains a specific tag ID."
  (cond ((atom expr) nil)
        ((and (eq (first expr) :tagged) (equal (second expr) id)) t)
        (t (some (lambda (x) (contains-tag-id-p x id)) expr))))

(defun apply-contraction-in-context (expr pair)
  "Locates the terms containing the pair's IDs and replaces them
   with the result of the physics rules."
  (let ((id1 (get-tag-id (first pair))) (id2 (get-tag-id (second pair))))
    (if (eq (first expr) '*)
        (let* ((args (cdr expr))
               (pos1 (position-if (lambda (x) (contains-tag-id-p x id1)) args))
               (pos2 (position-if (lambda (x) (contains-tag-id-p x id2)) args)))
          (if (= pos1 pos2)
              ;; Same Term (Recursive contraction inside one trace)
              (let ((new-args (copy-list args)))
                (setf (nth pos1 new-args) (apply-contraction-rules (nth pos1 args) pair))
                (cons '* new-args))
              ;; Different Terms (Merging two traces)
              (let* ((t1 (nth pos1 args)) (t2 (nth pos2 args))
                     (merged (apply-contraction-rules (list t1 t2) pair))
                     (others (loop for i from 0 for x in args unless (or (= i pos1) (= i pos2)) collect x)))
                (cons '* (cons merged others)))))
        ;; Base case: Single term
        (apply-contraction-rules expr pair))))

(defun apply-contraction-to-pairing (expr pairings)
  "Processes a specific pairing on the expression.
   1. Expands expression.
   2. Iteratively applies contractions.
   3. Handles branching sums generated by expansions."
  (let ((ex (expand-expression expr)))
    (cond
      ;; Case: Sum generated by expansion -> Map over all terms
      ((and (consp ex) (eq (first ex) '+))
       (cons '+ (loop for term in (cdr ex) 
                      for res = (apply-contraction-to-pairing term pairings)
                      if (and (consp res) (eq (first res) '+)) append (cdr res) else collect res)))
      ;; Case: No pairs left -> Return result
      ((null pairings) ex)
      ;; Case: Contract -> Apply next pair
      (t (apply-contraction-to-pairing (apply-contraction-in-context ex (first pairings)) (rest pairings))))))

;; ============================================================
;; 8. THE SAFE CLEANER
;; ============================================================

(defun simplify-expression (expr)
  "The 'Structural Sieve'.
   Does NOT perform complex algebra. It only cleans up the structure.
   
   Operations:
   1. Flattens nested Sums.
   2. Prunes dead branches in Products (* A 0 B -> 0).
   3. Prunes dead branches in Operators (Tr ... 0 ... -> 0).
   4. Applies safe atomic physics rules (e.g. Tr[Sigma3] -> 0)."
  (cond
    ((atom expr) expr)
    
    ;; === SUM (+) ===
    ((eq (first expr) '+)
     (let* ((args (mapcar #'simplify-expression (cdr expr)))
            (flat (loop for a in args if (and (consp a) (eq (first a) '+)) append (cdr a) else collect a))
            (clean (remove-if (lambda (x) (equal x 0)) flat)))
       (cond ((null clean) 0)
             ((= (length clean) 1) (first clean))
             (t (cons '+ clean)))))

    ;; === PRODUCT (*) ===
    ((eq (first expr) '*)
     (let* ((args (mapcar #'simplify-expression (cdr expr)))
            (flat (loop for a in args if (and (consp a) (eq (first a) '*)) append (cdr a) else collect a)))
       (cond
         ((member 0 flat) 0)
         (t (let ((clean (remove-if (lambda (x) (equal x 1)) flat)))
              (cond ((null clean) 1)
                    ((= (length clean) 1) (first clean))
                    (t (cons '* clean))))))))

    ;; === OPERATORS (TR / PARA / PERP) ===
    ((member (first expr) '(tr PARA PERP))
     (let ((args (mapcar #'simplify-expression (cdr expr))))
       (cond
         ((member 0 args) 0)
         ((and (eq (first expr) 'tr) (= (length args) 1) (eq (first args) 'SIGMA3)) 0)
         ((and (eq (first expr) 'PERP) (= (length args) 1) 
               (or (eq (first args) 'PI) (and (consp (first args)) (eq (first (first args)) 'PI)))) 0)
         (t (cons (first expr) args)))))

    (t (cons (first expr) (mapcar #'simplify-expression (cdr expr))))))

;; ============================================================
;; 9. MATHEMATICA INTERFACE
;; ============================================================
(defun to-mathematica-string (expr &optional context)
  "Recursively converts Lisp structures to Mathematica string syntax.
   Context: 'MATRIX (uses Dot . separator) or 'SCALAR (uses Times * separator)."
  (cond
    ((numberp expr) (format nil "~A" expr))
    ((symbolp expr)
     (case expr
       (t-coupling "t") (TR "Tr") (PARA "Para") (PERP "Perp") (ID "Id")
       (SIGMA3 "Subscript[\\[Sigma], 3]") (SIGMA1 "Subscript[\\[Sigma], 1]")
       (PI "\\[CapitalPi]") (PHI "\\[CapitalPhi]")
       (UBAREU "\\[CapitalEpsilon]") (VARPHI "\\[Phi]")
       ;; FIXED: Wrapped nil in a list ((nil) "") to satisfy CASE syntax
       ((nil) "") 
       (otherwise (string-capitalize (symbol-name expr)))))
    ((listp expr)
     (let ((op (first expr)) (args (rest expr)))
       (cond
         ((eq op :tagged) (format nil "Tagged[~A, ~A]"
                                  (get-tag-id expr) (to-mathematica-string (get-tag-content expr) 'MATRIX)))
         
         ;; Derivatives: (d x) -> Der[x]
         ((eq op 'd) (format nil "Der[~A]" (to-mathematica-string (first args) 'SCALAR)))
         
         ;; Arithmetic formatting
         ((eq op '+) 
          (if (null args) "0"
              (let ((clean (remove-if (lambda (s) (string= s "")) 
                                      (mapcar (lambda (x) (to-mathematica-string x context)) args))))
                (format nil "(~{~A~^ + ~})" clean))))
         
         ((eq op '-)
          (if (= (length args) 1) (format nil "(-~A)" (to-mathematica-string (first args) context))
              (format nil
                      "(~A - ~A)" (to-mathematica-string (first args) context) (to-mathematica-string (second args) context))))
         ((eq op '/)
          (if (= (length args) 1) (format nil "(1 / ~A)" (to-mathematica-string (first args) context))
              (format nil
                      "(~A / ~A)" (to-mathematica-string (first args) context) (to-mathematica-string (second args) context))))
         
         ;; Product (* vs .) based on Context
         ((eq op '*)
          (if (null args) (if (eq context 'MATRIX) "Id" "1")
              (let* ((sep (if (eq context 'MATRIX) " . " " * "))
                     (clean (remove-if (lambda (s) (string= s "")) 
                                       (mapcar (lambda (x) (to-mathematica-string x context)) args))))
                (format nil (concatenate 'string "(~{~A~^" sep "~})") clean))))
         
         ;; Power
         ((eq op 'expt) (format nil
                                "(~A^~A)" (to-mathematica-string (first args) context) (to-mathematica-string (second args) 'SCALAR)))

         ;; Propagator formatting with Dot Product for Derivatives
         ((eq op 'PI)
          (let ((coord (first args)) (derivs (rest args)))
            (if (null derivs) (format nil "\\[CapitalPi][~A]" (to-mathematica-string coord 'SCALAR))
                (format nil "(~{~A~^ . ~} . \\[CapitalPi][~A])" 
                        (mapcar (lambda (x) (to-mathematica-string x 'SCALAR)) derivs) (to-mathematica-string coord 'SCALAR)))))
         
         ;; Matrix Ops (Switch context to MATRIX)
         ((member op '(TR PARA PERP))
          (let ((clean (remove-if (lambda (s) (string= s "")) 
                                  (mapcar (lambda (x) (to-mathematica-string x 'MATRIX)) args))))
            (format nil "~A[~A]" (to-mathematica-string op) 
                    (format nil "~{~A~^ . ~}" clean))))
         
         ;; Standard functions (Keep SCALAR context)
         (t (format nil
                    "~A[~{~A~^, ~}]" (to-mathematica-string op) (mapcar (lambda (x) (to-mathematica-string x 'SCALAR)) args))))))))

(defun format-topology (pairing)
  "Helper to format the list of pairs into Mathematica syntax: {{1,2}, {3,4}}."
  (format nil "{~{~A~^, ~}}" 
          (loop for pair in pairing collect (format nil "{~D, ~D}" (second (first pair)) (second (second pair))))))

(defun print-simulation-object (sim-data)
  "Helper: Prints a single simulation object as a Mathematica Association."
  (let ((skeleton (getf sim-data :skeleton)) (diagrams (getf sim-data :diagrams)))
    (format t "<|~%  \"Skeleton\" -> ~A,~%  \"Diagrams\" -> {~%" (to-mathematica-string skeleton 'MATRIX))
    (loop for d in diagrams for i from 1 for is-last = (= i (length diagrams)) do
          (format t "    <| \"ID\" -> ~D, \"Status\" -> \"~A\", \"Topology\" -> ~A, \"Value\" -> ~A |>~A~%"
                  (getf d :id) (getf d :status) (format-topology (getf d :pairing))
                  (to-mathematica-string (getf d :result) 'SCALAR) (if is-last "" ",")))
    (format t "  }~%|>")))

(defun generate-mathematica-report (sim-data-list)
  "Iterates over the list of simulation results and prints a Mathematica List of Associations."
  (format t "{~%")
  (loop for sim in sim-data-list
        for i from 1
        for is-last = (= i (length sim-data-list)) do
        (print-simulation-object sim)
        (format t "~A~%" (if is-last "" ",")))
  (format t "}~%"))

;; ============================================================
;; 10. WICK CONTRACTION SOLVER
;; ============================================================

(defun collect-terms (expr)
  "Splits an expression tree into a flat list of additive terms.
   Input: (+ A B (+ C D)) -> (A B C D)
   Input: A -> (A)"
  (if (and (consp expr) (eq (first expr) '+))
      (cdr expr)
      (list expr)))

(defun solve-wick-contraction (expr field-names)
  (let* ((norm (normalize expr))
         (tagged (with-tagging (tag-fields norm field-names)))
         (expanded (simplify-expression (expand-expression tagged)))
         (terms (collect-terms expanded)))

    (loop for term in terms collect
          (let* ((items (extract-tags term))
                 (groups (group-tags-by-type items))
                 (n-fields (length items)))
            
            (list :skeleton term
                  :diagrams
                  (cond
                    ((zerop n-fields)
                     (list (list :id 1 :pairing nil :result term :status :survived)))
                    
                    ;; Check Parity of EACH group
                    ((some (lambda (entry) (oddp (length (cdr entry)))) groups)
                     (list (list :id 1 :pairing nil :result 0 :status :vanished)))
                    
                    (t (loop for p in (generate-stratified-pairing items) 
                             for i from 1
                             for res = (simplify-expression (apply-contraction-to-pairing term p))
                             collect (list :id i 
                                           :pairing p 
                                           :result res 
                                           :status (if (equal res 0) :vanished :survived))))))))))

;; ============================================================
;; 11. API INTEGRATION & CLI
;; ============================================================

(defun main (expr &key (rules *default-contraction-rules*) 
                       (fields *default-fields*))
  "The High-Level Entry Point.
   
   Arguments:
     expr: The S-expression to solve (e.g. '(Tr ...)).
     rules: (Optional) Alist mapping Field-Types to Contraction-Functions.
     fields: (Optional) List of Field-Types to tag and contract."
  
  (let ((*contraction-rules* rules))
    (generate-mathematica-report (solve-wick-contraction expr fields))))

(defun run-from-cli ()
  "Entry point for CLI execution (called by Mathematica).
   Reads arguments from stdin, evaluates, and prints the report."
  (let ((args (rest sb-ext:*posix-argv*)))
    (when args
      (let* ((expr (read-from-string (first args)))
             (result (eval expr)))
        (format t "~a" result)
        (force-output)
        (sb-ext:exit :code 0)))))

(run-from-cli)
