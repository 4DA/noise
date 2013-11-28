#lang racket
(require racket/gui)
(require math/matrix)
(require math/array)
(require racket/draw)
(require racket/flonum)

(struct cell (tl tr bl br) #:transparent)
(struct node (gr x y xy) #:transparent)
(struct grid (nodes xc yc) #:transparent)
(struct vec2 (x y) #:transparent)

(define (make-bmp-dc w h)
  (new bitmap-dc% [bitmap (make-object bitmap% w h)]))

(define (clamp val)
  (if (< val 0) 0
      (if (> val 255) (remainder (floor val) 255) (floor val))))

(define (clamp-byte b)
  (if (> b 255) 255 b))

(define (make-exp-noise-col int)
  (let ([gc (inexact->exact (floor (+ 127 (* int 127))))])
    (make-color (clamp (expt (* gc 4) 2)) (clamp (/ gc 3)) gc)))

(define (conv-to-0-255 fl)
  (inexact->exact (floor (+ 127 (* fl 127)))))

(define (make-gray-col int)
  (let ([gc (conv-to-0-255 int)])
    (make-color gc gc gc)))

(define (make-gray-argb int)
  (let ([gc (inexact->exact (abs  (floor (* int 127))))])
    (bytes 255 gc gc gc)))

(define (v-normalize v)
  (let* ([x (vec2-x v)] [y (vec2-y v)]
         [len (flsqrt (fl+ (fl* x x) (fl* y y)))])
    (vec2 (fl/ x len) (fl/ y len))))

(define (make-grid xsz ysz cell-side)
  (define (gen-urv3)
    (v-normalize 
     (vec2  (sub1 (* 2 (random))) (sub1  (* 2 (random))))))

  (let ([nodes-x (add1 (/ xsz cell-side))] 
        [nodes-y (add1 (/ ysz cell-side))])

    (grid (for*/fold ([cv (vector)]) ( [ny nodes-y] [nx nodes-x])
            (let ([xp (* nx cell-side)]
                  [yp (* ny cell-side)])
              (vector-append cv (vector 
                                 (node (gen-urv3) xp yp
                                       (vec2 (exact->inexact xp) (exact->inexact yp)))))))
     nodes-x
     nodes-y)))

(define (get-cell g num)
  (let* ([xsz (grid-xc g)] [ysz (grid-yc g)]
         [tl (+ num (floor (/ num (sub1 xsz))))] 
         [tr (add1 tl)] [bl (+ tl xsz)] 
         [br (add1 bl)] [nodes (grid-nodes g)])

    (cell (vector-ref nodes tl)
          (vector-ref nodes tr)
          (vector-ref nodes bl)
          (vector-ref nodes br))))

(define (render-noise bmp-dc xsz ysz cell-side)
  (let* ([pix-col (make-object color% 0 0 0)]
         [gr (make-grid xsz ysz cell-side)]
         [cx (sub1 (grid-xc gr))]
         [cy (sub1 (grid-yc gr))])

    (printf "Rendering ~ax~a cells...\n" cx cy)
    (for ([c (* cx cy)])
      (render-cell bmp-dc (get-cell gr c)))))

(define pix-col (make-object color% 0 0 0))

(define (lerp a b x)
  (fl+ (fl* a (fl- 1.0 x)) (fl* b x)))

(define (cos-interp a b x)
  (let* ([ft (fl* x 3.14159)]
         [f (fl* (fl- 1.0 (flcos ft)) 0.5)])
    (fl+ (fl* a (fl- 1.0 f)) (fl* b f))))

(define (q-poly t)
  (+ (fl* 6.0 (flexpt t 5.0))
     (- (fl* 15.0 (flexpt t 4.0)))
     (fl* 10.0 (flexpt t 3.0))))

(define (q-interp a b x)
  (let ([qpx (q-poly x)])
    (fl+ (fl* a qpx)
	 (fl* b (fl- 1.0 qpx)))))

;; u ---- v
;; |      |
;; |      |
;; s------t

(define (v-dot v1 v2)
  (fl+ (fl* (vec2-x v1) (vec2-x v2))
       (fl* (vec2-y v1) (vec2-y v2))))

(define (v- v1 v2)
  (vec2 (fl- (vec2-x v1) (vec2-x v2))
        (fl- (vec2-y v1) (vec2-y v2))))

(define cl-min 255)
(define cl-max 0)

(define (render-cell bmp cur-cell)
  (let ([cbl (cell-bl cur-cell)]
        [cbr (cell-br cur-cell)]
        [ctl (cell-tl cur-cell)]
        [ctr (cell-tr cur-cell)])

    (define (get-ifs x y)
      (define (get-inf cur-node norm-xy)

        (fl/ (v-dot (v- (vec2 x y)
			 norm-xy)
		     (node-gr cur-node))
	     1.0))
      
      (values (get-inf cbl (vec2 0.0 1.0)) 
              (get-inf cbr (vec2 1.0 1.0))
              (get-inf ctl (vec2 0.0 0.0)) 
              (get-inf ctr (vec2 1.0 0.0))))

    (let* ([xsi (node-x ctl)] [ysi (node-y ctl)]
           [xei (node-x ctr)] [yei (node-y cbl)]
           [cwi (- xei xsi)] [chi (- yei ysi)]
           [xs (->fl xsi)] [ys (->fl ysi)]
           [xe (->fl xei)] [ye (->fl yei)]
           [cw (->fl cwi)] [ch (->fl chi)])

      (for* ([yi (in-range ysi yei)]
             [xi (in-range xsi xei)])

        (let*-values ([(y) (->fl yi)]
                      [(x) (->fl xi)]
                      [(s t u v) 
                       (get-ifs (fl/ (fl- x xs) cw) (fl/ (fl- y ys) ch))])

          (let* ([xv (fl/ (- x xs) cw)] [yv (fl/ (- y ys) ch)]
                 [a (q-interp t s xv)] [b (q-interp v u xv)]
                 [res (q-interp a b yv)]
		 [colspc-pix (conv-to-0-255 res)])
	    (when (< colspc-pix cl-min) (set! cl-min colspc-pix))
	    (when (> colspc-pix cl-max) (set! cl-max colspc-pix))
            (send bmp set-pixel xi yi (make-gray-col res))))))
    (cons cl-min cl-max)))

(define (get-noise-octave buf1 f)
  (define xsz 512)
  (define ysz 512)

  (let* ([blen (* xsz ysz 4)]
         [buf2 (make-bytes blen)])

    ;; (send dc get-argb-pixels 0 0 xsz ysz buf1)
    ;; (send dc get-argb-pixels 0 0 xsz ysz buf2)

    (for* ([x xsz] [y ysz])
      (let* ([fx (remainder (* x f) xsz)] 
             [fy (remainder (* y f) ysz)]
             [soft (* 4 (+ fx (* fy ysz)))]  [doft (* 4 (+ x (* y ysz)))]
             [a (bytes-ref buf1 soft)]
             [r (bytes-ref buf1 (+ 1 soft))]
             [g (bytes-ref buf1 (+ 2 soft))]
             [b (bytes-ref buf1 (+ 3 soft))])
        
        (bytes-copy! buf2 doft (bytes a r g b))))
    buf2))

(define (merge-noise! mbuf buf-oct a)
  (define xsz 512)
  (define ysz 512)

  (define (mc m o)
    (inexact->exact (floor (+ (* m (- 1.0  a)) (* o a)))))

  (let* ([blen (* xsz ysz 4)]
	 [min 255]
	 [max 0])
    (for* ([x xsz] [y ysz])
      (let* ([oft (* 4 (+ x (* y ysz)))]
             [sr (bytes-ref mbuf (+ 1 oft))]
             [sg (bytes-ref mbuf (+ 2 oft))]
             [sb (bytes-ref mbuf (+ 3 oft))]
             [dr (bytes-ref buf-oct (+ 1 oft))]
             [dg (bytes-ref buf-oct (+ 2 oft))]
             [db (bytes-ref buf-oct (+ 3 oft))])
	
	(let ([rc (mc sr dr)]
	      [gc (mc sg dg)]
	      [bc (mc sb db)])


	  (bytes-copy! mbuf (+ 1 oft) 
		       (bytes rc
			      gc
			      bc)))))))

(define (gamma-correct buf min max)
  (define xsz 512)
  (define ysz 512)
  
  (let* ([blen (* xsz ysz 4)]
	 [rtm (/ 255 max)])
    (for* ([x xsz] [y ysz])
      (let* ([oft (* 4 (+ x (* y ysz)))]
             [sr (bytes-ref buf (+ 1 oft))]
             [sg (bytes-ref buf (+ 2 oft))]
             [sb (bytes-ref buf (+ 3 oft))])
	
	(bytes-copy! buf (+ 1 oft) 
		     (bytes (floor (* (- sr min) rtm))
			    (floor (* (- sg min) rtm))
			    (floor (* (- sb min) rtm))))))
    (cons min max)))

(define (run)
  (let* ([xsz 512] [ysz 512]
         [frame (new frame% [label "perlin noise"]
		     [width xsz]
		     [height ysz])]
         [src-frame (new frame% [label "perlin noise source"]
		     [width xsz]
		     [height ysz])]
	 [canvas (new canvas% [parent frame])]
	 [src-canvas (new canvas% [parent src-frame])]
	 [dc (send canvas get-dc)]
	 [src-dc (send src-canvas get-dc)]
	 [bmp-dc1 (make-bmp-dc xsz ysz)]
         [cell-side 64])

    (send frame show #t)
    (send src-frame show #t)

    (sleep/yield 1)

    (time
     (render-noise bmp-dc1 xsz ysz 
                   cell-side))
    (let* ([blen (* xsz ysz 4)]
           [mbuf (make-bytes blen)]
           [res-buf (make-bytes blen)])

      (send bmp-dc1 get-argb-pixels 0 0 xsz ysz mbuf)
      (send bmp-dc1 get-argb-pixels 0 0 xsz ysz res-buf)

      (gamma-correct mbuf cl-min cl-max)

      (merge-noise!  res-buf (get-noise-octave mbuf 2) 0.3)
      (merge-noise!  res-buf (get-noise-octave mbuf 4) 0.1)
      (merge-noise!  res-buf (get-noise-octave mbuf 8) 0.01)

      (send bmp-dc1 set-argb-pixels 0 0 xsz ysz res-buf)
      (send dc draw-bitmap (send bmp-dc1 get-bitmap) 0 0)

      (send bmp-dc1 set-argb-pixels 0 0 xsz ysz mbuf)
      (send src-dc draw-bitmap (send bmp-dc1 get-bitmap) 0 0))

    (printf "done.\n")
    frame))

(run)
