#lang racket

;;; Vadim Zaliva <vzaliva@cmu.edu>
;;; A. E. Eiben and J. E. Smith, Introduction to Evolutionary Computing (Natural Computing Series). Springer, 2007, p. 300.
;;; Excercise 3-7

(require srfi/1)

(define (random-string l)
  (let loop ((i l))
    (if (= i 0) '()
        (cons (random 2) (loop (sub1 i))))))

(define (random-strings n l)
  (let loop ((i n))
    (if (= i 0) '()
        (cons (random-string l) (loop (sub1 i))))))

(define (flip x) (if (= x 0) 1 0))

(define (mutate x)
  (if (null? x)
      '()
      (cons (if (< (random) PM)
                (flip (car x))
                (car x))
            (mutate (cdr x)))
      ))

;; "One Max" fitness function
(define (fitness x)  (reduce + 0 x))
;; Optimal solution check predicate
(define (optimal? i) (every (lambda (x) (= x 1)) i))

;; One-point crossover
(define (crossover a b)
  (let* ((r (add1 (random (- (length a) 2))))
         (ah (take a r)) (at (drop a r))
         (bh (take b r)) (bt (drop b r)))
    (list (append ah bt)
          (append bh at))
    ))

(define (recombine-and-mutate a b)
  (if (< (random) PC)
      (map mutate (crossover a b))
      (list a b)))

;; "roulette whell" algorithm
(define (roulette-wheel p)
  (let* ((pf (map fitness p))
         (tpf (reduce + 0 pf))
         (Psel (map (lambda (x) (/ x tpf)) pf))
         (a (let aloop ((al Psel) (v 0))
              (if (null? al)
                  '()
                  (let ((nv (+ v (car al))))
                    (cons nv (aloop (cdr al) nv)))
                  ))))
    (map (lambda (c)
           (let aselect ((al a) (pl p) (r (random)))
             (if (>= (car al) r)
                 (car pl)
                 (aselect (cdr al) (cdr pl) r))
             ))
         p)
    ))

;; Stochastic universal sampling algorithm
(define (SUS p)
  (let* ((pf (map fitness p))
         (tpf (reduce + 0 pf))
         (Psel (map (lambda (x) (/ x tpf)) pf))
         (a (let aloop ((al Psel) (v 0))
              (if (null? al)
                  '()
                  (let ((nv (+ v (car al))))
                    (cons nv (aloop (cdr al) nv)))
                  )))
         (mu (length p))
         (r0 (/ (random) mu)))
    (let loop ((r r0) (al a) (pl p))
      (if (null? pl)
          '()
          (if (<= r (car al))
              (cons (car pl) (loop (+ r (/ 1 mu)) al pl))
              (loop r (cdr al) (cdr pl)))))
    ))

(define (new-generation p)
  (let recombine-generation ((mp (SUS p)))
    (if (< (length mp) 2)
        mp
        (append (recombine-and-mutate (car mp) (cadr mp)) (recombine-generation (drop mp 2))))))

(define (evolve p0 max-gen)
  (let eloop ((g 0) (p p0))
    (if (= g max-gen)
        (values '() #f)
        (if (any optimal? p)
            (values '() #t)
            (let* ((pf (map fitness p))
                   (maxf  (exact->inexact (reduce max (car pf) pf)))
                   (minf  (exact->inexact (reduce min (car pf) pf)))
                   (avgf (exact->inexact (/ (reduce + 0 pf) (length pf)))))
              (let-values (((t f) (eloop (add1 g) (new-generation p))))
                (values (cons (list g minf avgf maxf) t) f)
                ))
            ))
    ))

;;;--- Model parametrization and execution ---

;; Initialize Random Number Generator to get repeatable results
(random-seed 1234)

;; Population size
(define L 25)
;; Crossover rate
(define PC 0.7)
;; Mutation rate
(define PM (/ 1 L))
;; Population size (mu)
(define PSIZE 100)
;; Max. number of generations
(define MAXGEN 10000)

;; Randomly initialize population
(define population0 (random-strings PSIZE L))

;; Do 10 runs
(let hloop ((a 1))
  (if (= a 11)
      #f
      (let-values (((t f) (evolve population0 MAXGEN)))
        (if f
            (printf "Run #~A: solution found in ~A generations\n" a (car (last t)))
            (printf "Run #~A: solution not found\n" a)
            )
        (hloop (add1 a))
        )))

;; Do an additional single run and record min, max, avg as CSV file for later plotting
(with-output-to-file "generations.csv"
  (lambda ()
    (begin
      (let-values (((t f) (evolve population0 MAXGEN)))
        (map (lambda (x) (printf "~A,~A,~A\n" (second x) (third x) (fourth x))) t))
      #t)))
