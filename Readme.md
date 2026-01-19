**A sample Mathematica notebook to calculate contractions!**

ClearAll["Global`*"];

SetDirectory[NotebookDirectory[]];

<< NLSMDriver`;

(*2. Define a Dummy Model*)

fullexpr = 
  1/(2 t) Tr[\[CapitalPhi][r] . Subscript[\[Sigma], 3] . W[r] . 
      W[r] . \[CapitalPhi][r] . Subscript[\[Sigma], 3] . W[r] . 
      W[r]] + \[Lambda]/(2 t)
     Tr[\[CapitalPhi][r] . Subscript[\[Sigma], 
      3] . \[CapitalPhi][r] . Subscript[\[Sigma], 3] . W[r] . W[r] . 
      W[r] . W[r]] + (1 + \[Lambda])/(2 t)
     Tr[\[CapitalPhi][r] . Subscript[\[Sigma], 3] . 
      W[r] . \[CapitalPhi][r] . Subscript[\[Sigma], 3] . W[r] . 
      W[r] . W[r]];

(*3. Run the Driver*)

NLSMDriver`GenerateDiagrams[fullexpr, "Report", "S2p", "Display" -> True]


**The model architecture**
```mermaid
graph TD
    User([User / Prototype Notebook])
    Driver[NLSMDriver<br/><i>The Manager</i>]
    Algebra[NLSMAlgebra<br/><i>Math Engine</i>]
    Lisp[NLSMLisp<br/><i>Lisp Interface & Translator</i>]
    Visuals[NLSMVisuals<br/><i>Pure Renderer</i>]
    ExtLisp((External<br/>Lisp Program: Contraction Engine))
    FileSystem[(File System)]

    User -->|1. Calls with full expression| Driver
    
    Driver -->|2. Splits Terms & Sends to LISP| Lisp
    Lisp <-->|CLI| ExtLisp
    
    Driver -->|3. Algebraically Simplifies Results| Algebra
    
    Driver -->|4. Sends Data for Diagrams & Report| Visuals
    Visuals -->|5. Exports PDF| FileSystem
    Visuals -.->|Return| User

    classDef main fill:#f9f,stroke:#333,stroke-width:2px;
    classDef sub fill:#e1f5fe,stroke:#333;
    class Driver main;
    class Algebra,Lisp,Visuals sub;
```

**---Module to be implemented---**

**Grand generator module**: This module will generate all interaction terms from a given action, restricted to a specific accuracy (e.g., 2-loop). After generating the full expansion, it will feed the final expression to the Driver. This module will essentially replace the current manual input role of the "User" for a complete automation pipeline.

**Integration module**: An added module between the Visuals and Driver. This will perform the necessary integration and keep the singular term (logarithmic quantum corrections in 2d). These will be fed back to the visuals for generating report.

Here is the current execution flow of the Program-->

```mermaid
sequenceDiagram
    participant U as User
    participant D as NLSMDriver
    participant L as NLSMLisp
    participant Li as External Lisp
    participant A as NLSMAlgebra
    participant V as NLSMVisuals

    U->>D: GenerateDiagrams(Expr)
    
    loop Every Term
        D->>D: IsolateFactor(Term)
        
        %% The Lisp Execution Chain
        D->>L: ToLispString(Trace-Expr)
        L->>Li: EvaluateModel()
        
        Note over Li: 1. Wick Contractions<br/>2. Identify Topology<br/>3. Compute Raw Value
        
        Li-->>L: Returns {Skeleton, Topologies, Value}
        L-->>D: Returns {Skeleton, Topologies, Value}
        
        %% The Math & Vis Chain
        D->>A: FinalSimplify(Factor * Value)
        A-->>D: Returns Simplified Expression
        
        D->>D: Format Title (TraditionalForm)
        
        D->>V: saveReport(EnhancedData)
        V-->>U: Saves PDF to Disk
    end
    
    D-->>U: Returns Visual Grid
```
