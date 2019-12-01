      subroutine imtql2(nm,n,d,e,z,ierr)

      integer i,j,k,l,m,n,ii,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag
      ierr = 0
      if (n .eq. 1) go to 1001

      do 100 i = 2, n
  100 e(i-1) = e(i)

      e(n) = 0.0e0

      do 240 l = 1, n
         j = 0
!     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = dabs(d(m)) + dabs(d(m+1))
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
  110    continue

  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         g = (d(l+1) - p) / (2.0e0 * e(l))
         r = pythag(g,1.0e0)
         g = d(m) - p + e(l) / (g + sign(r,g))
         s = 1.0e0
         c = 1.0e0
         p = 0.0e0
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            r = pythag(f,g)
            e(i+1) = r
            if (r .eq. 0.0e0) go to 210
            s = f / r
            c = g / r
            g = d(i+1) - p
            r = (d(i) - g) * s + 2.0e0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
!     .......... form vector ..........
            do 180 k = 1, n
               f = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * f
               z(k,i) = c * z(k,i) - s * f
  180       continue

  200    continue

         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0e0
         go to 105
!     .......... recover from underflow ..........
  210    d(i+1) = d(i+1) - p
         e(m) = 0.0e0
         go to 105
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)

         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue

         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p

         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
!
  300 continue

      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
