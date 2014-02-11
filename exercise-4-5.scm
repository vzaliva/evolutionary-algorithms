#lang racket

;;; Stochastic Optimization, Homework 2
;;; Vadim Zaliva <vzaliva@cmu.edu>

;;; A. E. Eiben and J. E. Smith, Introduction to Evolutionary Computing (Natural Computing Series). Springer, 2007, p. 300.
;;; Exercise 4-5

(require srfi/1)
(require (planet schematics/random:1:0/random))

(define ackley-min -30)
(define ackley-max 30)

(define (feasible-allele? a)
  (and (<= a ackley-max) (>= a ackley-min)))

;; Get variables from representation
(define (get-vars x) (take x N))
;; Get params from representation
(define (get-params x) (drop x N))

;; We use $x_i in [-30,30]$ input domain
(define (random-X)
  (- (* (random)  (- ackley-max ackley-min)) ackley-max))

;; Sigma initialization per:
;; T. Bäck and H. Schwefel, “An overview of evolutionary algorithms for parameter optimization,” Evol. Comput., no. 1, pp. 1–23, 1993.
(define (random-Sigma) 3)

(define (random-individual l)
  (append
   (let loopX ((i l))
     (if (= i 0) '()
         (cons (random-X) (loopX (sub1 i)))))
   (let loopS ((i l))
     (if (= i 0) '()
         (cons (random-Sigma) (loopS (sub1 i)))))
   ))

(define (feasible-individual? x)
  (every feasible-allele? x))

(define (random-population)
  (let loop ((i MU))
    (if (= i 0) '()
        (cons (random-individual N) (loop (sub1 i))))))

(define (simpleUnzip2 l)
  (let-values (((x s) (unzip2 l)))
    (append x s)))

;; Uncorrelated Mutation with n Step Sizes (sec. 4.4.2)
(define (mutate x)
  (let* ((xx (get-vars x))
         (xs (get-params x))
         (t0 (/ 1 (sqrt (* 2 (sqrt N)))))
         (t1 (/ 1 (sqrt (* 2 N))))
         (g1 (* t1 (random-gaussian))))
    (simpleUnzip2 (map (lambda (x0 s0)
                         (let* ((s1 (max E0 (* s0
                                               (exp (+ g1 (* t0 (random-gaussian)))))))
                                (x1 (+ x0 (* s1 (random-gaussian)))))
                           (list x1 s1)))
                       xx xs))))

(define (ackley  x0)
  (let ((x (take x0 N)))
    (+
     (* -20 (exp
             (* -0.2 (sqrt
                      (/ (reduce + 0 (map (lambda (i) (* i i)) x)) N)))))
     (* -1 (exp
            (/ (reduce + 0
                       (map (lambda (i) (cos (* 2 pi i))) x)
                       ) N)
            ))
     20
     (exp 1)
     )))

(define (random-pick l)
  (list-ref l (random (length l))))

;; Discrete recombination of 2 givel alleles
(define (discrete-recombination a b)
  (if (= (random 2) 0) a b))

;; Intermediate recombination of 2 givel alleles
(define (intermediate-recombination a b)
  (/ (+ a b) 2))

;; Local recombination strategy
(define (recombine-local fetcher combiner p)
  (map combiner
       (fetcher (random-pick p))
       (fetcher (random-pick p))))

;; Global recombination strategy
(define (recombine-global fetcher combiner p)
  (let rloop ((i 0))
    (if (= i N)
        '()
        (let ((p0 (fetcher (random-pick p)))
              (p1 (fetcher (random-pick p))))
          (cons
           (combiner (list-ref p0 i) (list-ref p1 i))
           (rloop (add1 i))))
        )))

;; discrete recombination of object variables and global intermediate
;; recombination of strategy parameters.
(define (recombine p)
  (append
   (recombine-local get-vars discrete-recombination p)
   (recombine-global get-params intermediate-recombination p)))

(define (new-child p)
  (let ((c (mutate (recombine p))))
    (if (feasible-individual? c)
        c
        (new-child p))))

(define (new-generation p)
  (take (sort
         (let loop ((i 0))
           (if (= i LA)
               '()
               (cons (new-child p) (loop (add1 i)))))
         fitness-order-p)
        MU))

(define (fitness-order-p p0 p1)
  (< (ackley p0) (ackley p1)))

(define (mean l) (/ (reduce + 0 l) (length l)))

;;;--- Model parametrization and execution ---

;; Initialize Random Number Generator to get repeatable results
(random-seed 1234)

(define MAXGEN 200000) ;; Max. number of generations
(define N 30) ;; Dimensionality of Ackley function
(define MU 30) ;; Population size
(define LA 200) ;; Selection pool size
(define E0 1E-10) ;; Min. bondary for sigma
(define MAXRUNS 100) ;; Number of runs

;; Termination criteria per:
;; T. Bäck and H. Schwefel, “An overview of evolutionary algorithms for parameter optimization,” Evol. Comput., no. 1, pp. 1–23, 1993.
(define ALMOST0 1E-4) ;; close enough to minimum

(define (run-once)
  (let gloop ((g 0)
              (p (random-population)))
    (let* ((pf (map ackley p))
           (f (fold min +inf.0 pf)))
      (printf "Generation #~A: Fb=~A, avg(F)=~A, max(|Xb|)=~A, avg(Sb)=~A\n" g f
              (mean pf)
              (fold max 0 (map abs (get-vars (car p))))
              (mean (get-params (car p))))
      (if (or (= g MAXGEN) (< f ALMOST0))
          (get-vars (car p))
          (gloop (add1 g) (new-generation p)))
      )))

(define (run-repeated n)
  (printf "*** Run #~A\n" n)
  (if (= n 0)
      '()
      (cons (run-once) (run-repeated (sub1 n)))))

;; Helper functions for statistics calculations
(define (v+ a b)  (map + a b))
(define (v- a b)  (map - a b))
(define (zerov n)  (make-list n  0))
(define (sq x) (* x x))

(define (main)
  (let* (
         (res (run-repeated MAXRUNS))
         (s (fold v+ (zerov N) res))
         (m (map (lambda (x) (/ x MAXRUNS)) s))
         (sd
          (map (lambda (z) (sqrt (/ z (- MAXRUNS 1))))
               (fold v+ (zerov N)
                     (map (lambda (x) (map sq (v- x m))) res)))
          ))
    (printf "Mean:\n~A\n" m)
    (printf "Std. dev:\n~A\n" sd)
    ))

(main)
