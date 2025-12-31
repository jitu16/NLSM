(defvar *counter* 0)

;; ============================================================
;; 1. NORMALIZATION & UTILS
;; ============================================================

(defun expand-trace-arg (term)
  "Expands an exponent expression into a list of repeated terms.
   Input: (expt X 3)
   Output: (X X X)
   If input is not an exponent, returns it as a single-item list."
  (if (and (consp term) (eq (first term) 'expt))
      (loop repeat (third term) collect (copy-tree (second term)))
      (list term)))

(defun normalize (expr)
  "Walks the entire expression tree to handle pre-processing tasks.
   Currently, it expands powers inside traces (Tr[X^2] -> Tr[X X])."
  (cond
    ((atom expr) expr)
    ;; Case: Trace -> Expand arguments using helper
    ((eq (first expr) 'tr)
     (cons 'tr (mapcan #'expand-trace-arg (cdr expr))))
    ;; Case: Generic List -> Recurse down
    (t (cons (first expr) (mapcar #'normalize (cdr expr))))))

(defun remove-nth (n list)
  "Removes the item at index N from the list (non-destructively)."
  (if (zerop n) (cdr list)
      (cons (car list) (remove-nth (1- n) (cdr list)))))

(defun deep-copy (expr) 
  "Creates a fresh deep copy of the expression tree."
  (copy-tree expr))

;; ============================================================
;; 2. TAGGING & PAIRING
;; ============================================================

(defun tag-fields (expr field-name)
  "Recursively walks the expression. Every time it finds the target field,
   it wraps it in a tagged structure with a unique ID.
   CRITICAL FIX: Uses (copy-tree expr) to ensure no shared memory references."
  (cond
    ;; Found the field -> Tag it and increment counter
    ((and (consp expr) (eq (first expr) field-name))
     (incf *counter*)
     (list :tagged *counter* (copy-tree expr))) ;; <--- COPIES STRUCTURE
    ;; Generic List -> Recurse
    ((consp expr)
     (cons (tag-fields (car expr) field-name)
           (tag-fields (cdr expr) field-name)))
    ;; Atom -> Return as is
    (t expr)))

(defmacro with-tagging (&body body)
  "Sets up the tagging environment (resets the counter)."
  `(let ((*counter* 0)) ,@body))

(defun extract-tags (expr)
  "Returns a flat list of all (:tagged ...) objects found in the tree.
   Used to build the list of items available for Wick contraction."
  (cond
    ((atom expr) nil)
    ((eq (first expr) :tagged) (list expr))
    (t (mapcan #'extract-tags expr))))

(defun generate-pairs (items)
  "Generates ALL possible Wick contraction pairings for the list of items.
   Returns a list of scenarios, where each scenario is a list of pairs.
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

;; ============================================================
;; 3. PHYSICS RULES (SPLIT & MERGE)
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
  "For Different-Trace contraction (Tr[A W]).
   Cyclically rotates the trace content so that W (id) is at the end,
   then returns everything else (A). Effectively computes A = Tr_content / W."
  (let ((pos (position-if (lambda (x) (and (consp x) (eq (first x) :tagged) (= (second x) id))) trace-content)))
    (append (subseq trace-content (1+ pos)) (subseq trace-content 0 pos))))

(defun get-splitting-coefficients (trace-content id1 id2)
  "For Same-Trace contraction (Tr[A W1 B W2]).
   Returns two values: the list A (between W1 and W2) and the list B (wrapping around).
   Crucial for the splitting formula."
  (let ((pos1 (position-if (lambda (x) (and (consp x) (eq (first x) :tagged) (= (second x) id1))) trace-content))
        (pos2 (position-if (lambda (x) (and (consp x) (eq (first x) :tagged) (= (second x) id2))) trace-content)))
    (if (< pos1 pos2)
        (values (subseq trace-content (1+ pos1) pos2)
                (append (subseq trace-content (1+ pos2)) (subseq trace-content 0 pos1)))
        (values (append (subseq trace-content (1+ pos1)) (subseq trace-content 0 pos2))
                (subseq trace-content (1+ pos2) pos1)))))

(defun apply-physics-rules (skeleton pairing)
  "The Core Physics Logic.
   1. Extracts derivatives from the paired fields (W r dx).
   2. Constructs the Propagator (Pi ...).
   3. Determines if we are Splitting (same trace) or Merging (different traces).
   4. Returns the mathematical expression for the contraction result.
   FIXED: Deep copies extracted lists to prevent shared references."
  (let* ((tag1 (first pairing)) (tag2 (second pairing))
         (id1 (second tag1)) (id2 (second tag2))
         (content1 (third tag1)) (content2 (third tag2))
         
         ;; --- Extract Derivatives ---
         ;; Input format expected: (W r (d x) (d y)...)
         (r1 (second content1)) (r2 (second content2))
         ;; Use cddr to skip 'W' and 'r', collecting all (d ...) terms
         (all-derivs (append (cddr content1) (cddr content2)))
         
         ;; --- Construct Propagator ---
         (propagator (list* 'PI (list '- r1 r2) all-derivs))
         
         ;; --- Locate Traces ---
         (trace1 (find-trace-containing skeleton id1))
         (trace2 (find-trace-containing skeleton id2)))

    (cond
      ((or (null trace1) (null trace2)) (error "Tags missing in skeleton."))
      
      ;; Case A: SPLITTING (Same Trace object)
      ;; Formula: -t/4 * (Tr[A_para Pi]Tr[B_para] - Tr[A_para S3 Pi]Tr[B_para S3])
      ((eq trace1 trace2)
       (multiple-value-bind (A-list B-list) 
           (get-splitting-coefficients (get-content-without-tr trace1) id1 id2)
         
         ;; COPY THE LISTS HERE to ensure independence
         (let ((A-para (apply-op 'PARA (deep-copy A-list)))
               (B-para (apply-op 'PARA (deep-copy B-list))))
           `(* (/ (- t-coupling) 4)
               (- (* (tr ,@A-para ,propagator) (tr ,@B-para))
                  (* (tr ,@(deep-copy A-para) sigma3 ,(deep-copy propagator)) 
                     (tr ,@(deep-copy B-para) sigma3)))))))

      ;; Case B: MERGING (Different Trace objects)
      ;; Formula: -t/2 * Tr[A_perp Pi B_perp]
      (t 
       (let* ((A-list (get-single-coefficient (get-content-without-tr trace1) id1))
              (B-list (get-single-coefficient (get-content-without-tr trace2) id2))
              
              ;; COPY THE LISTS HERE to ensure independence
              (A-perp (apply-op 'PERP (deep-copy A-list)))
              (B-perp (apply-op 'PERP (deep-copy B-list))))
         `(* (/ (- t-coupling) 2) (tr ,@A-perp ,propagator ,@B-perp)))))))

;; ============================================================
;; 4. EXPANSION ENGINE (THE ARCHITECT)
;; ============================================================

(defun has-tagged-field-p (expr)
  "Checks if an expression tree contains any tagged fields.
   Used to decide if we need to expand an operator to free a field."
  (cond ((atom expr) nil)
        ((eq (first expr) :tagged) t)
        (t (some #'has-tagged-field-p expr))))

(defun distribute-linear-op (op args)
  "Distributes a linear operator (Tr, Para, Perp) over a Sum.
   Also handles flattening products inside traces if necessary.
   Input: (Tr (+ A B)) -> (+ (Tr A) (Tr B))"
  (let ((pos-of-sum (position-if (lambda (x) (and (consp x) (eq (first x) '+))) args)))
    (if (null pos-of-sum) (cons op args)
        (let ((sum-terms (cdr (nth pos-of-sum args)))
              (prefix (subseq args 0 pos-of-sum))
              (suffix (subseq args (1+ pos-of-sum))))
          (cons '+ (loop for term in sum-terms collect
                         (let* ((middle (if (and (eq op 'tr) (consp term) (eq (first term) '*))
                                            (cdr term) (list term))))
                           (expand-expression (cons op (append (deep-copy prefix) middle (deep-copy suffix)))))))))))

;; --- HELPER: PREVENTS EMPTY OPERATORS ---
(defun make-op-term (op content)
  "Creates (OP content...). 
   CRITICAL FIX: 
     If content is empty:
     PARA -> 'ID
     PERP -> 0 (This kills invalid terms like Para[W])"
  (if (null content) 
      (if (eq op 'PARA) 'ID 0)
      (cons op content)))

(defun expand-operator (op terms)
  "Expands (OP A W B) based on the parity of W (assumed Odd/Perp).
   Recursive expansion ensures multiple Ws are handled correctly.
   Deep copies ensure no shared memory."
  (if (not (has-tagged-field-p terms)) 
      (cons op (deep-copy terms)) ;; Base case: No fields, just copy.
      
      (let* ((pos (position-if (lambda (x) (and (consp x) (eq (first x) :tagged))) terms))
             (W (nth pos terms))
             (L (subseq terms 0 pos))
             (R (subseq terms (1+ pos))))

        ;; W is inherently PERP (Odd).
        ;; Multiplication Rules:
        ;; PARA (Even) -> Requires L * W * R to be Even.
        ;;    Since W is Odd, L*R must be Odd.
        ;;    So L and R must have DIFFERENT parities (Para/Perp or Perp/Para).
        ;;
        ;; PERP (Odd)  -> Requires L * W * R to be Odd.
        ;;    Since W is Odd, L*R must be Even.
        ;;    So L and R must have SAME parities (Para/Para or Perp/Perp).

        (if (eq op 'PARA)
            ;; Parent is PARA: Split into Mixed Parity
            `(+ (* ,(make-op-term 'PARA (deep-copy L)) 
                   ,(deep-copy W) 
                   ,(expand-expression (make-op-term 'PERP (deep-copy R)))) ;; Recurse on R
                (* ,(make-op-term 'PERP (deep-copy L)) 
                   ,(deep-copy W) 
                   ,(expand-expression (make-op-term 'PARA (deep-copy R)))))

            ;; Parent is PERP: Split into Same Parity
            `(+ (* ,(make-op-term 'PARA (deep-copy L)) 
                   ,(deep-copy W) 
                   ,(expand-expression (make-op-term 'PARA (deep-copy R)))) ;; Recurse on R
                (* ,(make-op-term 'PERP (deep-copy L)) 
                   ,(deep-copy W) 
                   ,(expand-expression (make-op-term 'PERP (deep-copy R)))))))))


(defun expand-expression (expr)
  "The Main Expansion Routine.
   Ensures that tagged fields are brought to the top level of products.
   Handles:
   - Converting (-) to (* -1)
   - Flattening Sums
   - Distributing operators over sums
   - Expanding Para/Perp using the structural rules."
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
;; 5. ORCHESTRATOR
;; ============================================================

(defun get-tag-id (item) (second item))

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
                (setf (nth pos1 new-args) (apply-physics-rules (nth pos1 args) pair))
                (cons '* new-args))
              ;; Different Terms (Merging two traces)
              (let* ((t1 (nth pos1 args)) (t2 (nth pos2 args))
                     (merged (apply-physics-rules (list t1 t2) pair))
                     (others (loop for i from 0 for x in args unless (or (= i pos1) (= i pos2)) collect x)))
                (cons '* (cons merged others)))))
        ;; Base case: Single term
        (apply-physics-rules expr pair))))

(defun process-scenario (expr pairings)
  "Processes a specific pairing scenario on the expression.
   1. Expands expression to make fields accessible.
   2. Iteratively applies contractions from the pairings list.
   3. Handles branching sums generated by expansions."
  (let ((ex (expand-expression expr)))
    (cond
      ;; Case: Sum generated by expansion -> Map over all terms
      ((and (consp ex) (eq (first ex) '+))
       (cons '+ (loop for term in (cdr ex) 
                      for res = (process-scenario term pairings)
                      if (and (consp res) (eq (first res) '+)) append (cdr res) else collect res)))
      ;; Case: No pairs left -> Return result
      ((null pairings) ex)
      ;; Case: Contract -> Apply next pair
      (t (process-scenario (apply-contraction-in-context ex (first pairings)) (rest pairings))))))

;; ============================================================
;; 6. THE SAFE CLEANER (The Sieve)
;; ============================================================

(defun simplify-expression (expr)
  "The 'Structural Sieve'.
   Does NOT perform complex algebra. It only cleans up the structure to prevent
   zero-diagrams from being passed to Mathematica.
   
   Operations:
   1. Flattens nested Sums: (+ A (+ B)) -> (+ A B).
   2. Prunes dead branches in Products: (* A 0 B) -> 0.
   3. Prunes dead branches in Operators: (TR A 0 B) -> 0.
   4. Applies safe atomic physics rules: Tr[Sigma3] -> 0, Perp[Pi] -> 0."
  (cond
    ((atom expr) expr)
    
    ;; === SUM (+) ===
    ;; Flattens nested sums and removes zeros.
    ((eq (first expr) '+)
     (let* ((args (mapcar #'simplify-expression (cdr expr)))
            (flat (loop for a in args if (and (consp a) (eq (first a) '+)) append (cdr a) else collect a))
            (clean (remove-if (lambda (x) (equal x 0)) flat)))
       (cond ((null clean) 0)
             ((= (length clean) 1) (first clean))
             (t (cons '+ clean)))))

    ;; === PRODUCT (*) ===
    ;; Flattens nested products, removes identity elements (1), and annihilates if 0 is present.
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
    ;; Applies structural cleanup and safe atomic physics rules.
    ((member (first expr) '(tr PARA PERP))
     (let ((args (mapcar #'simplify-expression (cdr expr))))
       (cond
         ;; If ANY argument is 0, the linear operator evaluates to 0.
         ((member 0 args) 0)
         
         ;; Atomic Physics Rule: Tr[Sigma3] vanishes.
         ((and (eq (first expr) 'tr) (= (length args) 1) (eq (first args) 'SIGMA3)) 0)
         
         ;; Atomic Physics Rule: Perp[Pi] vanishes because Pi is diagonal.
         ((and (eq (first expr) 'PERP) (= (length args) 1) 
               (or (eq (first args) 'PI) (and (consp (first args)) (eq (first (first args)) 'PI)))) 0)
         
         (t (cons (first expr) args)))))

    ;; === GENERIC LIST ===
    (t (cons (first expr) (mapcar #'simplify-expression (cdr expr))))))



;; ============================================================
;; 7. MATHEMATICA INTERFACE
;; ============================================================

(defun to-mma-string (expr &optional context)
  "Recursively converts Lisp structures to Mathematica string syntax.
   Context: 'MATRIX (uses Dot . separator) or 'SCALAR (uses Times * separator).
   UPDATED: Strictly filters empty strings to prevent ' . . ' syntax errors."
  (cond
    ((numberp expr) (format nil "~A" expr))
    ((symbolp expr)
     (case expr
       (t-coupling "t") (TR "Tr") (PARA "Para") (PERP "Perp") (ID "Id")
       (SIGMA3 "Subscript[\\[Sigma], 3]") (PI "\\[CapitalPi]") (PHI "\\[CapitalPhi]")
       (UbarEU "\\[CapitalEpsilon]") (varphi "\\[Phi]")
       (nil "") ;; Handle NIL by returning empty string
       (otherwise (string-capitalize (symbol-name expr)))))
    ((listp expr)
     (let ((op (first expr)) (args (rest expr)))
       (cond
         ;; Tagged items -> print value inside
         ((eq op :tagged) (format nil "Tagged[~A, ~A]" (first args) (to-mma-string (second args) 'MATRIX)))
         
         ;; Derivatives: (d x) -> Der[x]
         ((eq op 'd) (format nil "Der[~A]" (to-mma-string (first args) 'SCALAR)))
         
         ;; Arithmetic Formatting
         ((eq op '+) 
          (if (null args) "0"
              (let ((clean (remove-if (lambda (s) (string= s "")) 
                                      (mapcar (lambda (x) (to-mma-string x context)) args))))
                (format nil "(~{~A~^ + ~})" clean))))
         
         ((eq op '-) (if (= (length args) 1) (format nil "(-~A)" (to-mma-string (first args) context))
                         (format nil "(~A - ~A)" (to-mma-string (first args) context) (to-mma-string (second args) context))))
         ((eq op '/) (if (= (length args) 1) (format nil "(1 / ~A)" (to-mma-string (first args) context))
                         (format nil "(~A / ~A)" (to-mma-string (first args) context) (to-mma-string (second args) context))))
         
         ;; Product (* vs .) based on Context
         ((eq op '*)
          (if (null args) (if (eq context 'MATRIX) "Id" "1")
              (let* ((sep (if (eq context 'MATRIX) " . " " * "))
                     (clean (remove-if (lambda (s) (string= s "")) 
                                       (mapcar (lambda (x) (to-mma-string x context)) args))))
                (format nil (concatenate 'string "(~{~A~^" sep "~})") clean))))
         
         ;; Power
         ((eq op 'expt) (format nil "(~A^~A)" (to-mma-string (first args) context) (to-mma-string (second args) 'SCALAR)))

         ;; Propagator formatting with Dot Product for Derivatives
         ((eq op 'PI)
          (let ((coord (first args)) (derivs (rest args)))
            (if (null derivs) (format nil "\\[CapitalPi][~A]" (to-mma-string coord 'SCALAR))
                (format nil "(~{~A~^ . ~} . \\[CapitalPi][~A])" 
                        (mapcar (lambda (x) (to-mma-string x 'SCALAR)) derivs) (to-mma-string coord 'SCALAR)))))
         
         ;; Matrix Ops (Switch context to MATRIX)
         ((member op '(TR PARA PERP))
          (let ((clean (remove-if (lambda (s) (string= s "")) 
                                  (mapcar (lambda (x) (to-mma-string x 'MATRIX)) args))))
            (format nil "~A[~A]" (to-mma-string op) 
                    (format nil "~{~A~^ . ~}" clean))))
         
         ;; Standard functions (Keep SCALAR context)
         (t (format nil "~A[~{~A~^, ~}]" (to-mma-string op) (mapcar (lambda (x) (to-mma-string x 'SCALAR)) args))))))))

(defun format-topology-mma (pairing)
  "Helper to format the list of pairs into Mathematica syntax: {{1,2}, {3,4}}."
  (format nil "{~{~A~^, ~}}" 
          (loop for pair in pairing collect (format nil "{~D, ~D}" (second (first pair)) (second (second pair))))))

(defun generate-mma-report (sim-data)
  "Generates the final Mathematica Association string containing the Skeleton and all Diagrams."
  (let ((skeleton (getf sim-data :skeleton)) (diagrams (getf sim-data :diagrams)))
    (format t "<|~%  \"Skeleton\" -> ~A,~%  \"Diagrams\" -> {~%" (to-mma-string skeleton 'MATRIX))
    (loop for d in diagrams for i from 1 for is-last = (= i (length diagrams)) do
          (format t "    <| \"ID\" -> ~D, \"Status\" -> \"~A\", \"Topology\" -> ~A, \"Value\" -> ~A |>~A~%"
                  (getf d :id) (getf d :status) (format-topology-mma (getf d :pairing))
                  (to-mma-string (getf d :result) 'SCALAR) (if is-last "" ",")))
    (format t "  }~%|>~%")))

;; ============================================================
;; 8. MAIN & CLI
;; ============================================================

(defun solve-wick-contraction (expr field-name)
  "Top-level function to solve the model.
   1. Normalizes input.
   2. Tags fields.
   3. Generates pairs.
   4. Processes every scenario through expansion, contraction, and cleaning.
   5. Returns the structured simulation object."
  (let* ((norm (normalize expr))
         (tagged (with-tagging (tag-fields norm field-name)))
         (items (extract-tags tagged)))
    (if (oddp (length items)) (error "Odd number of fields.")
        (list :skeleton tagged :diagrams
              (loop for p in (generate-pairs items) for i from 1
                    for res = (simplify-expression (process-scenario tagged p))
                    collect (list :id i :pairing p :result res :status (if (equal res 0) :vanished :survived)))))))

(defun main (raw-input) 
  "Main entry point for local testing."
  (generate-mma-report (solve-wick-contraction raw-input 'W)))

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



