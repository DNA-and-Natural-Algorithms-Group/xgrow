
;; by Erik Winfree, July 1999
;;    added xgrow format output, July 2002

;; one can go through this file, using ESC-control-x in Emacs-Lisp mode,
;; to execute the commands one by one.  Be sure to have two windows open;
;; the second one will be used as space for the simulation

;; a tile is '(name north-bd east-bd south-bd west-bd)
;; there is a global function that determines binding domain strength from
;; from the bd name

(defun tilestrength (bd) "Tile System: binding domain strength"
  (length bd)
)
;; how many characters are the tile names?
(setq tilenamelen 1)

(defun tilespace (X Y) "Tile System: Init"
 (let ((x 0) (y 0) (m tilenamelen))
  (goto-char 0)
  (insert " ")
  (setq x 0) (while (< x (* X m)) (setq x (+ x 1)) (insert "-") ) 
  (insert " \n") 
  (setq y 0) 
  (while (< y Y) 
    (setq y (+ y 1)) 
    (insert "|")
    (setq x 0) (while (< x (* X m)) (setq x (+ x 1)) (insert " ") ) 
    (insert "|\n") 
  )
  (insert " ")
  (setq x 0) (while (< x (* X m)) (setq x (+ x 1)) (insert "-") ) 
  (insert " \n") 
 )
)

;; s is a tile (c N E S W) where c is a tilenamelen-character string

(defun tileplace (X Y x y s) "Tile System: Placement"
  (goto-char (+ (* y (+ (* X tilenamelen) 3)) 2 (* (- x 1) tilenamelen)))
  (delete-char tilenamelen) (insert (car s))
)

(defun string-after (pos len) 
 "Return length-LEN string in current buffer starting with character
  at position POS.  If POS+LEN-1 is out of range, the value is nil."
  (if (= len 0) 
      ""
      (if (null (char-after pos))
          nil
          (concat (char-to-string (char-after pos)) 
                  (string-after (+ pos 1) (- len 1)) 
          )
      )
  )
)

(string-after 1034 6)

;; tiles is a list of tiles; find one whose char matches location

(defun tileread (X Y x y tiles) "Tile System: Who's there?"
  (let* ( (pos (+ (* y (+ (* X tilenamelen) 3)) 2 (* (- x 1) tilenamelen)))
          (ts tiles) 
          (s (string-after pos tilenamelen)) )
    (while (not (or (null ts) (equal (car (car ts)) s) ))
      (setq ts (cdr ts)) )
    (if (null ts) nil (car ts)))
)

;; a tile is (c N E S W): char, binding domain symbols 
;; a tile set is a list of tiles
;; sites lists is e.g. ( (tile1 x1 y1) (tile2 x2 y2) )
;; site is (x1 y1)

(defun tilezap (sites site) "Tile System: remove location"
 (apply 'append
  (mapcar (lambda (s) (if (equal (cdr s) site) nil (list s))) sites)
 )
) 

(tilezap '( (a 1 2) (b 2 3) (c 4 5) ) '(2 3) )

(setq ABtiles '( ("A" "bb" "bb" "aa" "aa") ("B" "aa" "aa" "bb" "bb") ) )


;; debugging tool -- the action is always in the other window
(defun other-do (arg)
  (other-window 1)
  (print (eval arg))
  (other-window 1)
)

(other-do '(tilespace 40 10) )
(other-do '(tileplace 40 10 20 5 '("A" blah)) )
(other-do '(tileread 40 10 20 5 '(("G" blah)) ) )
(other-do '(tileread 40 10 20 5 ABtiles ) )

;; note: strength of binding domain symbols = length as strings
;; choose all tiles that could 2-bind in this site
(defun tilecompat (tiles bN bE bS bW) 
     "Tile System: select if bonds to neighbors?"
 (if (null tiles) tiles
  (let ( (sN (if (equal bN (car (cdr (car tiles)))) 
                 (tilestrength bN) 0))
         (sE (if (equal bE (car (cdr (cdr (car tiles))))) 
                 (tilestrength bE) 0))
         (sS (if (equal bS (car (cdr (cdr (cdr (car tiles)))))) 
                 (tilestrength bS) 0))
         (sW (if (equal bW (car (cdr (cdr (cdr (cdr (car tiles))))))) 
                 (tilestrength bW) 0))
       )
    (if (> (+ sN sE sS sW) 1)
        (cons (car tiles) (tilecompat (cdr tiles) bN bE bS bW))
        (tilecompat (cdr tiles) bN bE bS bW)
    )
  )
 )
)

;; should be faster
(defun tilecompat (tiles bN bE bS bW) 
     "Tile System: select if bonds to neighbors?"
 (if (null tiles) tiles
  (apply 'append (mapcar
   (lambda (s) 
     (let ( (sN (if (equal bN (car (cdr s))) 
                    (tilestrength bN) 0))
            (sE (if (equal bE (car (cdr (cdr s)))) 
                    (tilestrength bE) 0))
            (sS (if (equal bS (car (cdr (cdr (cdr s))))) 
                    (tilestrength bS) 0))
            (sW (if (equal bW (car (cdr (cdr (cdr (cdr s))))))
                    (tilestrength bW) 0))
          )
          (if (> (+ sN sE sS sW) 1) (list s) nil)
     )
   )
   tiles
  ))
 )
)

(tilecompat ABtiles "aa" nil nil nil)
(tilecompat ABtiles  nil nil "bb" nil)

;; which tiles can fit at site x y?  must be empty
(defun tilefit (tiles X Y x y) "Tile System: find possible tiles"
 (let ((tC (tileread X Y x y tiles))
       (bN (car (cdr (cdr (cdr (tileread X Y x (- y 1) tiles))))))
       (bE (car (cdr (cdr (cdr (cdr (tileread X Y (+ x 1) y tiles)))))))
       (bS (car (cdr (tileread X Y x (+ y 1) tiles))))
       (bW (car (cdr (cdr (tileread X Y (- x 1) y tiles))))) )
   (if (null tC) (tilecompat tiles bN bE bS bW) nil) 
 )
)


(other-do '(tileread 40 10 20 5 ABtiles) )
(other-do '(tilefit ABtiles 40 10 20 4) )

;; sites lists is e.g. ( (tile1 x1 y1) (tile2 x2 y2) )
;; site is (x1 y1)

(defun tileupdate (X Y x y tiles sites) 
     "Tile System: re-figure who can go here"
  (let ( (st (tilezap sites (list x y)))
         (ts (tilefit tiles X Y x y)) )
    (while (not (null ts))
       (setq st (cons (list (car ts) x y) st))
       (setq ts (cdr ts))
    )
    st
  )
)

(defun tileglue (X Y x y tile tiles sites) "Tile System: put tile & update"
 (let ((st sites))
  (tileplace X Y x y tile)
  (setq st (tileupdate X Y x y tiles st))
  (if (< x X) (setq st (tileupdate X Y (+ x 1) y tiles st)))
  (if (< y Y) (setq st (tileupdate X Y x (+ y 1) tiles st)))
  (if (> x 1) (setq st (tileupdate X Y (- x 1) y tiles st)))
  (if (> y 1) (setq st (tileupdate X Y x (- y 1) tiles st)))
  st
 )
)

;; final sim: choose random from sites.  tileglue.  repeat.

(defun tilechoose (sites) "Tile System: should be random, but isn't"
  (car sites)
)

(defun tilegrow (X Y x y tiles tileseed) "Tile System: grow from seed"
  (let ((sites nil))
    (other-window 1) 
    (tilespace X Y) 
    (setq sites (tileglue X Y x y tileseed tiles sites))
    (while (not (null sites))
      (read-string (concat (prin1-to-string sites) " -- ?"))
      (let ((site (tilechoose sites)))
        (setq sites (tileglue X Y 
           (car (cdr site)) (car (cdr (cdr site))) (car site)
           tiles sites))
      )
    )
    (other-window 1)
  )
)

(car ABtiles)
(tilegrow 10 10 5 5 ABtiles (car ABtiles))


(setq CounterTiles 
  '( ("S" "rr" "-" "-" "ll")
     ("0" "0" "n" "0" "n") 
     ("1" "1" "n" "1" "n")
     ("O" "0" "c" "1" "c")
     ("l" "1" "c" "0" "n")
     ("R" "rr" "-" "rr" "c")
     ("L" "0" "ll" "-" "ll")    ) )

(tilegrow 20 20 19 19 CounterTiles (car CounterTiles))

(setq CombTiles
 '( ("1" "-" "12" "aa" "<")
    ("2" "-" "23" "aa" "12")
    ("3" "-" "34" "aa" "23")
    ("4" "-" ">"  "aa" "34")
    ("A" "aa" "<" "bb" ">")
    ("B" "bb" "<" "cc" ">")
    ("C" "cc" "<" "_" ">") ))

(tilegrow 20 20 10 10 CombTiles (car CombTiles))

(setq UnaryTiles
 '( ("1" "-" "12" "o"  "<")
    ("2" "-" "23" "o"  "12")
    ("3" "-" "34" "o"  "23")
    ("4" "-" "45" "aa" "34")
    ("5" "-" ">"  "o"  "45")
    ("O" "o" "o" "o" "o")
    ("A" "aa" "o" "o" "b")
    ("B" "o" "b" "aa" "o") ))

(tilegrow 20 20 10 10 UnaryTiles (car UnaryTiles))


(setq UnaryTiles2
 '( ("1" "-" "12" "r"  "<")
    ("2" "-" "23" "r"  "12")
    ("3" "-" "34" "r"  "23")
    ("4" "-" "45" "aa" "34")
    ("5" "-" ">"  "o"  "45")
    ("A" "aa" "o" "o" "b")
    ("B" "r" "b" "aa" "r")
    ("O" "o" "o" "o" "o")
    ("o" "r" "r" "r" "r") ))


(tilegrow 20 20 10 10 UnaryTiles2 (car UnaryTiles2))



(setq BinaryTiles
 '( ("#" "L" "12" "p"  "<")
    ("$" "0" "23" "p"  "12")
    ("%" "1" "34" "p"  "23")
    ("^" "0" "45" "AA" "34")
    ("&" "ZZ" "o" "o"  "45")
    ("A" "AA" "o" "o" "b")
    ("B" "p" "b" "AA" "p")
    ("O" "o" "o" "o" "o")
    ("o" "p" "p" "p" "p") 
    ("a" "b" "o" "o" "aa")
    ("b" "p" "aa" "b" "p")
    ("0" "0" "x" "0" "x")  
    ("1" "1" "x" "1" "x")
    ("c" "0" "c" "1" "c")
    ("C" "1" "c" "0" "n")
    ("N" "1" "n" "1" "n")
    ("n" "0" "n" "0" "n")
    ("r" "rr" "p" "r" "x")
    ("R" "RR" "p" "R" "x")
    ("K" "r" "aa" "ZZ" "c")
    ("s" "r" "p" "RR" "c")
    ("S" "R" "p" "rr" "n")
    ("l" "l" "x" "ll" "-")
    ("L" "L" "x" "LL" "-")
    ("t" "FC" "c" "L" "-")
    ("T" "LL" "c" "l" "-")
    ("u" "ll" "n" "l" "-")
    ("U" "LL" "n" "L" "-") ))

(tilegrow 50 50 5 40 BinaryTiles (car BinaryTiles))

;; modify global function so we can use 3 and 4 char bd names
(defun tilestrength (bd) (- (length bd) 2))

;; modify global variable so we can use two-character tile names
(setq tilenamelen 3)

(setq TMtiles1
 '( ("SS " "S__N" "S__E" "S__S" "S__W")

    ("NW " "0_N"  "  N"  "  W"  "0_W")
    ("NE " "0_N"  "0_E"  "  E"  "  N")
    ("SE " "  E"  "0_E"  "0_S"  "  S")
    ("SW " "  W"  "  S"  "0_S"  "0_W")

    ("IN " "A0_N" "  N"  "S__N" "  N")
    ("IE " "  E"  "A0_E" "  E"  "S__E")
    ("IS " "S__S" "  S"  "A0_S" "  S")
    ("IW " "  W"  "S__W" "  W"  "A0_W") ))

(tilegrow 20 20 10 10 TMtiles1 (car TMtiles1))

(other-do '(tileread 20 20 9 9 TMtiles1) )
(other-do '(tileplace 20 20 7 7 (car TMtiles1)) )

;; now add the TM rules, automatically making 4 variants

(defun tilemod (tile s) "Tile System: rename tile & bds w/ postfix s"
  (let ( (n (car tile))
         (N (car (cdr tile)))
         (E (car (cdr (cdr tile))))
         (S (car (cdr (cdr (cdr tile)))))
         (W (car (cdr (cdr (cdr (cdr tile)))))) )
    (list (concat n s) (concat N s) (concat E s) (concat S s) (concat W s))
  )
)

(defun tilerotate (tile) "Tile System: make 4 rotations of tile"
  (let ( (n (car tile))
         (N (car (cdr tile)))
         (E (car (cdr (cdr tile))))
         (S (car (cdr (cdr (cdr tile)))))
         (W (car (cdr (cdr (cdr (cdr tile)))))) )
    (list (tilemod (list n N E S W) "N")
          (tilemod (list n W N E S) "E")
          (tilemod (list n S W N E) "S")
          (tilemod (list n E S W N) "W")
   )
  )
)

(setq TMprototiles
 '( ("A0" "1_"  "-B" "A0_" "  ") ;; in state A, read 0; write 1, go left to B
    ("A1" "1_"  "  " "A1_" "-C")
    ("B0" "1_"  "  " "B0_" "-A")
    ("B1" "1_"  "-B" "B1_" "  ")
    ("C0" "1_"  "  " "C0_" "-B")
    ("C1" "***" "  " "C1_" "  ") 

    ("a0" "A0_" "  " "0_"  "-A") ;; arriving from left in state A; read 0
    ("0a" "A0_" "-A" "0_"  "  ")
    ("a1" "A1_" "  " "1_"  "-A")
    ("1a" "A1_" "-A" "1_"  "  ")
    ("b0" "B0_" "  " "0_"  "-B")
    ("0b" "B0_" "-B" "0_"  "  ")
    ("b1" "B1_" "  " "1_"  "-B")
    ("1b" "B1_" "-B" "1_"  "  ")
    ("c0" "C0_" "  " "0_"  "-C")
    ("0c" "C0_" "-C" "0_"  "  ")
    ("c1" "C1_" "  " "1_"  "-C")
    ("1c" "C1_" "-C" "1_"  "  ")

    (" 0" "0_"  "  " "0_"  "  ")  ;; maintain tape
    (" 1" "1_"  "  " "1_"  "  ")
  )
)

(setq joe (apply 'append (mapcar 'tilerotate TMprototiles) ))
(setq joe (cdr joe))

(setq TMtiles2 (apply 'append (mapcar 'tilerotate TMprototiles) ))
(setq TMtiles (append TMtiles1 TMtiles2))

(tilegrow 40 40 20 20 TMtiles (car TMtiles))

(length TMtiles)
(setq max-lisp-eval-depth 1000)

;; add code to write out the tile set in xgrow format to other window?





