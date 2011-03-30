C
C 6_0607  commented time-dependent output
C
      MODULE for77test

      contains

      SUBROUTINE t92(mode, a0,c0,row, nvar,neq,nonz, b, x,
     1               d,rr,v,vt,cna0,barm,
     2          mxp,ipd,stab,bar,stpbar,epsin,epsout,maxit,iter,
     3          a,nac,c,r,nrend,h,diag,kstep, merr)

      INCLUDE 'impl.inc'
      INTEGER*4 mode(4),nvar,neq,mxp,ipd,maxit,iter,merr(3),row(neq+1),h(9,-1:neq),kstep
     0INTEGER*4 c0(nonz),c(0:nac),r(0:nrend)
      REAL*8    stab,bar,stpbar,epsin,epsout
      REAL*8    vt(neq)
       DIMENSION a0(nonz),b(neq),x(nvar),d(nvar),rr(neq),v(neq)
       DIMENSION cna0(neq),barm(neq),a(nac),diag(nvar)
*                  РЕШЕНИЕ РАЗРЕЖЕННОЙ СИСТЕМЫ

*IO     mode   - выбор режима:
*IO              mode(1) = 0  - новая серия
*I               mode(2) = 0  - производится построчное неявное масштаб.
*I               mode(2) = 1  - производится упрощенное масштабирование
*I               mode(3) < 0  - печатается инф. о расходящейся попытке
*I               mode(3) = 0  - не печатается сборная информация
*I               mode(3) > 0  - не печатается сообщение об ошибках
*I               mode(4) = 0  - уравнение решается
*                               методом простых итераций
*I               mode(4) = 1  - решается сопряженное уравнение
*                               методом простых итераций
*I               mode(4) = 11 - уравнение решается
*                               методом бисопряженных градиентов
*I               mode(4) = 21 - решается сопряженное уравнение
*                               методом бисопряженных градиентов
*I      a0     - массив коэффициентов исходной системы,
*                расположенных построчно
*I      c0     - список столбцовых индексов коэффициентов из a0
*I      row    - список начал строк коэффициентов из a0
*I      nvar   - верхняя грань индексов неизвестных
*I      neq    - верхняя грань индексов уравнений
*I      nonz   - количество ненулевых коэффициентов в a0
*I      b      - правая часть системы уравнений
* O     x      - получаемое решение системы уравнений
* O     d      - рабочий вектор-добавка к решению
* O     rr     - рабочий вектор-невязка
* O     v      - рабочий вектор-правая часть нормированной системы
* O     vt     - рабочий вектор-для накопления
*       cna0   - рабочий вектор-C-нормы строк
*       barm   - рабочий вектор-барьеры
*I      mxp    - кол-во строк-кандидатов для выбора ведущего элемента
*I      ipd    - процент просветов.       0  =< ipd  =< 100
*IO     stab   - допуск устойчивости.    0.0 =< stab =< 1.0
*IO     bar    - барьер отсечения по 0
*IO     stpbar - шаг поиска оптимального барьера bar
*I      epsin  - требуемая  относительная точность в равномерной норме
* O     epsout - получаемая относительная точность в равномерной норме
*I      maxit  - наибольшее разрешенное количество итераций
* O     iter   - потребовавшееся количество итераций
*.O     a      - результат разложения исходной матрицы a0
*I      nac    - размер массивов a и c
*.O     c      - список столбцовых индексов коэффициентов из a
*.O     r      - список строчных индексов коэффициентов из a
*I      nrend  - размер массива r
*.O     h      - массив, задающий структуру разложения
*.O     diag   - массив ведущих элементов разложения
*.O     kstep  - число шагов декомпозиции ( = числу ведущих элементов )
* O     merr   - массив ошибок и статистических сведений
*                merr(1) = 0 - система решена. В этом случае:
*                              merr(2) = количеству сборок мусора в  c,
*                              merr(3) = количеству сборок мусора в  r.
*                merr(1) = 1 - не хватает памяти в  c и a
*                merr(1) = 2 - не хватает памяти в  r
*                merr(1) = 3 - даже вычиcления с машинной точностью
*                              не дают решения c заданной точностью.
*                              Необходимы вычисления с большей
*                              разрядной сеткой.
*                merr(1) =-1 - nzmax < nonz
*                merr(1) =-2 - nvar =< 0
*                merr(1) =-3 - neq  =< 0 или neq # nvar, если выбран
*                                     метод бисопряженных градиентов
*                merr(1) =-4 - nonz =< 0
*                merr(1) =-5 - nrend < nonz + 3
*                merr(1) =-6 - nac   < nonz
*                merr(1) =-7 - в  c0  есть индекс вне ( 1..nvar )
*                              а именно: это индекс  c(merr(2))
*                merr(1) =-8 - в  r0  есть индекс вне ( 1..neq )
*                              а именно: это индекс  r(merr(2))
                             

      REAL      t(3)
c,tmrast
c!     -        ,abig,one,z
      INTEGER*4 num,memor,ktry,i,jvar,ieq
      SAVE      num,memor
      
      DATA      num,memor / 0 , 0 /

*            правим неграмотно заданные параметры
      IF(stab.GT.1.0) stab=1.0/stab
      IF(stab.LT.0.0) stab=0.25
      IF(stpbar.LE.0.0) stpbar=1.0
      IF(stpbar.LT.1.0) stpbar=1.0/stpbar

      IF(mode(4).GE.11 .AND. mode(4).LE.29) THEN
         IF(nvar.NE.neq) THEN
            merr(1)=-3
            GOTO 90
         ENDIF
      ENDIF

      jvar=nvar
      ieq=neq
      nnz=row(ieq+1)-row(1)

      one=1.
      t(1)=tmrast(0)
      ktry=0
      num=num+1

      DO i=1,neq
        v(i)=b(i)
      END DO
      t(2)=tmrast(0)
      t(3)=t(2)
*                          вычисление C-норм строк
      DO i=1,neq
         abig=0.
         DO j=row(i),row(i+1)-1
            abig=max(abs(a0(j)),abig)
         ENDDO
         cna0(i)=abig
      ENDDO
      IF(mode(2).NE.0) THEN  ! вычисление C-нормы матрицы a0
         abig=0.
         DO i=1,neq
            abig=max(cna0(i),abig)
         ENDDO
      ENDIF
      IF(mode(1).NE.0) GOTO 30
      
*         размещение матрицы a0,
*         хранящейся в построчном формате,
*         в массиве  a  с целью дальнейшей декомпозиции

 20   CALL t90st4(a0,nnz,c0,row,nvar,neq,
     1           a,nac,c,r,nrend,h,ipd,merr)
      IF(merr(1).NE.0) GOTO 90

*            декомпозиция
                               
*                   попытка увеличить барьер
      IF(ktry.EQ.0) bar=bar*stpbar**.2
      IF(bar.GT.1.) bar=1.
      ktry=ktry+1
      IF(mode(2).EQ.0) THEN      ! определение барьеров
         DO i=1,neq
            barm(i)=cna0(i)*bar
         ENDDO
      ELSE                       ! упрощенные барьеры
         z=abig*bar
         DO i=1,neq
            barm(i)=z
         ENDDO
      ENDIF
      CALL t87d(a,nac,c,r,nrend,nvar,neq,row(neq+1)-row(1),
     1          h,diag,barm,stab,mxp,ipd,kstep,merr)
     

      t(2)=tmrast(0)
      t(3)=t(2)
      CALL mmr(h,nvar,neq,memor)
      IF(merr(1).NE.0) GOTO 90
      
*      если  матрица выродилась, попытаемся уменьшить барьер
*              ( не реализовано )

*            итерационное уточнение

 30   IF(mode(4).GE.11 .AND. mode(4).LE.25) THEN
         l0=0
         IF(mode(4).GE.21) l0=1
         itol=mod(mode(4),10)
         IF(itol.LE.0 .OR. itol.GT.4) itol=1
         CALL linbcg(a0,c0,row, nvar,neq,row(neq+1)-row(1), nvar,
     1               l0,itol,epsin,maxit,b, x,iter,epsout,
     2               a,nac,c,kstep,h,diag)
         IF(epsin.LT.epsout) iter=0
      ELSEIF(mode(4).EQ.0) THEN
         CALL t87ii6(nvar,neq,v,x,rr,d, vt,cna0,barm,
     1            maxit,epsin,epsout,iter,
     2            a,nac,c,kstep,h,diag,  a0,row(neq+1)-row(1),c0,row)
      ELSE
         CALL t87iit(nvar,neq,v,x,rr,d, vt,cna0,barm,
     1            maxit,epsin,epsout,iter,
     2            a,nac,c,kstep,h,diag,  a0,row(neq+1)-row(1),c0,row)
      ENDIF
     
      t(3)=tmrast(0)
      IF(iter.GT.0) GOTO 80
      
*            итерации расходятся, попытка уменьшения барьера
                                                           
 40   IF(mode(3).LT.0) CALL pr88
     -  (2,t,memor,merr,iter,nvar,neq,kstep,epsout,bar,stpbar,num,ktry,
     -nzmax,c0)
      IF(ktry.EQ.0) GOTO 20
      IF(mod(ktry,3).EQ.0) stpbar=2*stpbar
      IF(128.*bar+1.0 .EQ. one) GOTO 92
      bar=bar/stpbar
      GOTO 20

*                       посчитано

 80   mode(1)=1
      IF(mode(3).NE.0) CALL pr88
     -  (1,t,memor,merr,iter,nvar,neq,kstep,epsout,bar,stpbar,num,ktry,
     -nzmax,c0)
      IF(ktry.EQ.0) RETURN
      IF(stpbar.GE.8.) stpbar=stpbar/2
      RETURN
      
*                      ошибки

 90   IF(mode(3).GT.0) RETURN
      CALL pr88
     -  (3,t,memor,merr,iter,nvar,neq,kstep,epsout,bar,stpbar,num,ktry,
     -nzmax,c0)
      RETURN
 92   merr(1)=3
      IF(mode(3).GT.0) RETURN
      CALL pr88
     -  (4,t,memor,merr,iter,nvar,neq,kstep,epsout,bar,stpbar,num,ktry,
     -nzmax,c0)
      RETURN
      END SUBROUTINE

      SUBROUTINE t87d(a,nac,c,r,nrend,nvar,neq,nonz,
     1                h,diag,bar,stab,mxp,ipd,kstep,merr)
     
*        декомпозиция разреженной матрицы:    A = L * D * U
*        с отсечением малых элементов по масшт. массиву баръеров

      INCLUDE 'impl.inc'
      PARAMETER (maxp=100)
      INTEGER*4 nvar,neq,mxp,ipd,kstep
      INTEGER*4 nac,nrend,nonz,merr(3),h(9,-1:nvar),
     -          lk,mk,li,mi,lj,mj,
     -          ll,kz,rend,rend2,cf,rf,dc,dr,m2,mar,p,jc,l
     0INTEGER*4 c(0:nac),r(0:nrend),
     -          ind(2,maxp),
     -          i,j,s,t,m,n,jl,nk,ip,mp,ki,ig,ibig,ir,lin,is,js, iz2
*!      REAL*8    stab
*!     -          pd,barjer,zero
      DIMENSION a(nac),diag(nvar),   bar(neq),!
     -          ab(maxp)
*!     -          ,b,oa,aij,w
      LOGICAL   ek,ei

      zero=0.
      iz2=0
      m=neq
      n=nvar
      ll=nac+1
      kz=nonz
      rend=nrend-1
*!      barjer=bar
      ibig=m
      IF(n.GT.m) ibig=n
      mp=mxp
      IF(mp.LT.1) mp=1
      IF(mp.GT.maxp) mp=maxp
      pd=ipd*.01
      IF(pd.LT.0.) pd=0.
      IF(pd.GT.1.) pd=1.
      m2=2147483647 ! =2**31-1 ! =ibig*ibig
      rend2=rend-1
      merr(1)=0
      merr(2)=0
      merr(3)=0
      c(0)=1
      cf=h(2,m)+1
      rf=h(4,n)+1
      r(rend+1)=0
      r(rend)=1
      r(0)=1
      dc=pd*(ll-kz)/m
      dr=pd*(rend-kz)/n
      DO i=1,m
         h(8,i)=0
      ENDDO
      DO j=1,n
         h(7,j)=0
         h(9,j)=0
         diag(j)=0.
      ENDDO
 
*   Заполнение упорядоченной структуры строк

      h(7,-1)=0
      h(7,0)=0
      jl=0
      DO i=1,m                         ! включение в структуру
         ig=h(2,i)-h(1,i)              ! i-ой строки
         IF(h(7,ig).NE.0) THEN
            s=h(7,ig)
            h(6,i)=s
            j=h(5,s)
            h(5,i)=j
            h(6,j)=i
            h(5,s)=i
         ELSE
            h(7,ig)=i
            h(5,i)=i
            h(6,i)=i
         ENDIF
      ENDDO
                                       

*   200-499   Основной внешний цикл

      kstep=0
      nk=n
      
*   200-299   Выбор ведущего элемента (is,js), js=c(jc)

 200  IF(jl.NE.iz2) jl=jl-1            ! установка стрелки jl
 205  IF(h(7,jl).EQ.0) THEN            ! на начало наименьшей строки
         IF(jl.EQ.nk) GOTO 500
         jl=jl+1
         GOTO 205
      ENDIF
      kstep=kstep+1
      nk=n-kstep
      
*   по 260. Вычисление наименьшей допустимой
*           цены марковица mar

      mar=m2
      i=0
      s=0
      ki=0
      t=jl-1
      DO 260 ip=1,mp
         IF(i.NE.s) GOTO 220           ! следующая строка i
 210     IF(t.EQ.nk) GOTO 269
         t=t+1
         i=h(7,t)
         IF(i.EQ.0) GOTO 210
         s=i
 220     b=0.                          ! найдем b= равномерной
         li=h(1,i)                     ! норме строки
         mi=h(2,i)
         DO l=li,mi
            IF(ABS(b).LT.ABS(a(l))) b=a(l)
         ENDDO
         b=ABS(b)
         oa=b*stab                     ! oa=границе устойчивости
*            в строке вычисляем наименьшую длину
*            допустимого столбца
         ig=ibig
         DO l=li,mi
            IF(ABS(a(l)).GE.ABS(oa)) THEN
               j=c(l)
               lin=h(4,j)-h(3,j)
               IF(lin.LT.ig) ig=lin
            ENDIF
         ENDDO
         IF(ig*t-mar) 250,255,260      ! сравнение со старой ценой
 250        mar=ig*t                   ! начать новый отсчет
            ki=0
 255        ki=ki+1                    ! добавить в список
            ab(ki)=b                   ! строк-кандидатов
            ind(1,ki)=i
            ind(2,ki)=ig
 260        i=h(6,i)                   ! пропустить


*   находим максимальный /строчно-нормированный/
*   элемент с ценой Марковица  mar
 269     oa=-1.
         DO ip=1,ki
            i=ind(1,ip)
            ig=ind(2,ip)
            b=0.
*              из элементов строки с ценой mar выбираем наибольший
            DO l=h(1,i),h(2,i)
               j=c(l)
               IF(h(4,j)-h(3,j).EQ.ig .AND. ABS(a(l)).GE.ABS(b)) THEN
                  b=a(l)
                  p=l
               ENDIF
            ENDDO
*              сравниваем с предыдущими строками
            b=ABS(b)    ! b=ABS(b/ab(ip)) - если строчно-нормированный
            IF(b.GT.oa) THEN
               oa=b
               is=i
               jc=p
            ENDIF
         ENDDO
 280  CONTINUE
      js=c(jc)
      
*      CALL tpr(1,kstep,nvar,neq,is,js,a,nac,c,r,rend,h,diag)
      
      lk=h(1,is)
      lin=h(2,is)-lk
      IF(oa.NE.zero) oa=-1./a(jc)          ! выборка ведущего элемента
      a(jc)=a(lk)
      c(jc)=c(lk)
      c(lk)=js
 300  IF(cf.LT.ll-lin-1) GOTO 302          ! надо ли чистить c ?
      dc=pd*(cf-kz)/(m-kstep+2)            !      да
      CALL grbc(a,cf,c,m,h,dc,merr)
      IF(merr(1).EQ.1) GOTO 900
      GOTO 300
*         удаление из столбцов элементов is-ой строки
*         и перепись is-ой строки в конец c
 302  lk=h(1,is)
      mk=h(2,is)
      p=ll-mk-1
      DO jc=lk,mk
         j=c(jc)
         c(jc+p)=j
         c(jc)=0
         a(jc+p)=a(jc)
         diag(j)=a(jc)
         mj=h(4,j)
         l=h(3,j)
 310     IF(r(l).NE.is) THEN
            l=l+1
            GOTO 310
         ENDIF
         r(l)=r(mj)
         r(mj)=0
         h(4,j)=mj-1
      ENDDO
      diag(js)=-oa
      kz=kz-lin-1
      i=h(5,is)                        ! удаление is-ой строки из
      IF(i.EQ.is) THEN                 ! упорядоченной структуры строк
         h(7,lin)=0
      ELSE
         j=h(6,is)
         h(6,i)=j
         h(5,j)=i
         IF(h(7,lin).EQ.is) h(7,lin)=j
      ENDIF
      lk=ll-lin                         ! переопределение атрибутов
      h(1,is)=lk                        ! is-ой строки
      mk=ll-1
      h(2,is)=mk
      c(lk-1)=0
      h(8,is)=kstep
      ek=lk.EQ.ll
      IF(oa.EQ.zero) THEN                ! пропускаем строку,
         ll=lk                           ! если она нулевая
         h(8,is)=-1                      ! (но не пустая)
         GOTO 205
      ENDIF


*   по 499.   элементарные преобразования строк

      lin=h(4,js)-h(3,js)
      ll=lk-lin-1
      IF(lin.EQ.-1) GOTO 489
      DO 470 ir=0,lin
         l=h(3,js)+ir
         i=r(l)
         ig=h(2,i)-h(1,i)
 410     IF(cf+dc+nk.LT.ll) GOTO 412     ! нужна ли сборка мусора в c ?
         dc=pd*(cf-kz)/(m-kstep+2)       !          да
         CALL grbc(a,cf,c,m,h,dc,merr)
         IF(merr(1).EQ.1) GOTO 900
         GOTO 410
 412     s=h(5,i)                        ! удаление i-ой строки из
         IF(s.EQ.i) THEN                 ! упорядоченной структуры
            h(7,ig)=0                    ! строк
         ELSE
            j=h(6,i)
            h(6,s)=j
            h(5,j)=s
            IF(h(7,ig).EQ.i) h(7,ig)=j
         ENDIF
         li=h(1,i)                        ! найдем элемент a(i,js)
         jc=li
         mi=h(2,i)
 420     IF(c(jc).NE.js) THEN
            jc=jc+1
            GOTO 420
         ENDIF
         w=oa*a(jc)                       ! перепись и удаление
         a(ll+ir)=w                       ! элемента (i,js)
         c(ll+ir)=i
         c(jc)=c(mi)
         a(jc)=a(mi)
         c(mi)=0
         ei=li.EQ.mi
         mi=mi-1
         IF(ek) GOTO 460
         IF(ei) GOTO 440
         
*                проход по строке i
         DO jc=li,mi
            j=c(jc)
            IF(diag(j).NE.zero) THEN
               a(jc)=a(jc)+w*diag(j)
               diag(j)=0.
            ENDIF
         ENDDO
 440     CONTINUE
 
*                   проход по строке is

         barjer=bar(i)
         DO 450 jc=lk,mk
            j=c(jc)
            IF(diag(j).EQ.zero) THEN
               diag(j)=a(jc)
            ELSE
               aij=w*diag(j)
*                      сравнение с барьером
               IF(ABS(aij).LE.barjer) GOTO 450
               
*                      запись нового элемента в список строк
               IF(c(mi+1).NE.iz2) GOTO 610       ! есть ли место справа?
                  mi=mi+1                        !    да
                  a(mi)=aij
                  c(mi)=j
                  GOTO 650
 610           IF(c(li-1).NE.iz2) GOTO 620       ! есть ли место слева?
                  li=li-1                        !    да
                  a(li)=aij
                  c(li)=j
                  GOTO 650
 620           cf=cf+dc
               p=cf-li                           ! записать строку i
               DO l=li,mi                        !     в конец
                  a(l+p)=a(l)
                  c(l+p)=c(l)
                  c(l)=0
               ENDDO
               li=cf
               mi=mi+p+1
               a(mi)=aij
               c(mi)=j
               cf=mi+1
               
*                      запись нового элемента в список столбцов
 650           mj=h(4,j)
               IF(r(mj+1).NE.iz2) GOTO 660       ! есть ли место внизу?
                 mj=mj+1                         !    да
                 r(mj)=i
                 h(4,j)=mj
                 GOTO 450
 660           lj=h(3,j)                         ! есть ли место вверху?
               IF(r(lj-1).NE.iz2) GOTO 670
                 lj=lj-1                         !    да
                 r(lj)=i
                 h(3,j)=lj
                 GOTO 450
 670           IF(r(rf).NE.iz2) THEN             ! подправим
                  rf=rf+1                        ! "свободно" r
                  GOTO 670
               ENDIF
                  
*                      переписать столбец j в конец списка r
 680           rf=rf+dr
               p=rf-lj
               IF(mj+p.LT.rend2) GOTO 682  ! хватит ли места для копии?
               dr=pd*(rend-kz)/(nk+2)      !          нет
               CALL grbr(rend,rf,r,n,h,dr,merr)
               IF(merr(1).EQ.2) GOTO 900
               p=rf
               lj=h(3,j)
               mj=h(4,j)
               GOTO 680
 682           IF(lj.GT.mj) GOTO 692       !          да
               DO l=lj,mj
                  r(l+p)=r(l)
                  r(l)=0
               ENDDO
 692           mj=mj+p+1
               r(mj)=i
               h(4,j)=mj
               h(3,j)=rf
               rf=mj+1
            ENDIF
 450     CONTINUE
         IF(mi.GE.cf) cf=mi+1
 460     h(1,i)=li
         h(2,i)=mi
         kz=kz-ig
         ig=mi-li
         kz=kz+ig
         
*   включение i-ой строки в упорядоченную структуру строк
         IF(h(7,ig).ne.0) THEN
            s=h(7,ig)
            h(6,i)=s
            j=h(5,s)
            h(5,i)=j
            h(6,j)=i
            h(5,s)=i
         ELSE
            h(7,ig)=i
            h(5,i)=i
            h(6,i)=i
         ENDIF
C        h(7,ig)=i     ! если включаем в начало
 470  CONTINUE
      
*   подготовка diag к следующему шагу

      lj=h(3,js)
      mj=h(4,js)
      DO l=lj,mj
         r(l)=0
      ENDDO
 489  IF(ek) GOTO 499
         DO jc=lk,mk
c ! убрать след. строку, если отказываемся от нормировки верхнего тр-ка
            a(jc)=oa*a(jc)                  ! нормировка верхнего тр-ка
            diag(c(jc))=0.
         ENDDO
 499  h(9,js)=kstep
      h(3,js)=ll
      h(4,js)=lk-1
      GOTO 200
      
      
      
*   обратная перестановка

 500  CONTINUE
 
*      CALL tpr(2,kstep,nvar,neq,is,js,a,nac,c,r,rend,h,diag)
 
      DO i=1,m
         h(7,i)=h(8,i)
      ENDDO
      DO i=1,m
         j=h(7,i)
         IF(j.GT.0) h(8,j)=i
      ENDDO
      DO j=1,n
         h(7,j)=h(9,j)
      ENDDO
      DO j=1,n
         i=h(7,j)
         IF(i.GT.0) h(9,i)=j
      ENDDO
      RETURN
   
*               не хватает c или r

 900  RETURN
      END SUBROUTINE
      

      SUBROUTINE grbc(a,cf,c,m,h,dc,merr)
      
*             сборка мусора в c

      INCLUDE 'impl.inc'
     0INTEGER*4 m
      INTEGER*4 cf,dc,h(9,-1:m),merr(3),
     -          cf1,dc1,l,p,q
     0INTEGER*4 c(0:cf),
     -          s,t
      DIMENSION a(cf)
      
      merr(2)=merr(2)+1
      cf1=cf-1
 805  dc1=dc-1
      DO 810 s=1,m                        ! проставляем метки в начале
         IF(h(8,s).NE.0) GOTO 810         ! каждой оставшейся непустой
         q=h(1,s)                         ! строки
         t=h(2,s)-q
         IF(t.EQ.-1) GOTO 810
         h(2,s)=t
         h(1,s)=c(q)
         c(q)=-s
 810  CONTINUE
      p=1                                 ! просмотр c
      DO 828 q=1,cf1
         IF(c(q)) 820,828,825
 820        IF(q-p.LE.dc) GOTO 822        ! начало строки
            DO l=p,p+dc1                  ! вставим dc нулей
               c(l)=0
            ENDDO
            p=p+dc
 822        s=-c(q)
            c(q)=h(1,s)
            h(1,s)=p
            h(2,s)=p+h(2,s)
 825        c(p)=c(q)                     ! перепись q --> p
            a(p)=a(q)
            p=p+1
 828  CONTINUE
      IF(p.EQ.cf) GOTO 840                ! уплотнились ли?
         DO q=p,cf1                       ! да
            c(q)=0                        ! чистим хвост
         ENDDO
         cf=p
      RETURN
 840  IF(dc.EQ.0) GOTO 910                ! уплотнения нет
      dc=0
      GOTO 805
                      
*               не хватает c

 910  merr(1)=1
      END SUBROUTINE
      

      SUBROUTINE grbr(rend,rf,r,n,h,dr,merr)
      
*             сборка мусора в r
      
     0INTEGER*4 n
      INTEGER*4 rend,rf,dr,h(9,-1:n),merr(3),
     -          rf1,dr1,l,p,q
     0INTEGER*4 r(0:rend),
     -          s,t
     
      merr(3)=merr(3)+1
      IF(rf.GT.rend) rf=rend
      rf1=rf-1
 855  dr1=dr-1
      DO 860 s=1,n                      ! проставляем метки в начале
         IF(h(9,s).NE.0) GOTO 860       ! каждого оставшегося непустого
         q=h(3,s)                       ! столбца
         t=h(4,s)-q
         IF(t.EQ.-1) GOTO 860
         h(4,s)=t
         h(3,s)=r(q)
         r(q)=-s
 860  CONTINUE
      p=1                               ! просмотр r
      DO 878 q=1,rf1
         IF(r(q)) 870,878,875
 870        IF(q-p.LE.dr) GOTO 872      ! начало столбца
            DO l=p,p+dr1                ! вставим dr нулей
               r(l)=0
            ENDDO
            p=p+dr
 872        s=-r(q)
            r(q)=h(3,s)
            h(3,s)=p
            h(4,s)=p+h(4,s)
 875        r(p)=r(q)                   ! перепись q --> p
            p=p+1
 878  CONTINUE
      IF(p.EQ.rf) GOTO 890              ! уплотнились ли?
         DO q=p,rf1                     ! да
            r(q)=0                      ! чистим хвост
         ENDDO
         rf=p
      RETURN
 890  IF(dr.EQ.0) GOTO 920              ! уплотнения нет
      dr=0
      GOTO 855
      
*               не хватает r

 920  merr(1)=2
      END SUBROUTINE
      

      SUBROUTINE t87ii6(nvar,neq,b,x,r,d, t,cna0,rm,
     1                 maxit,epsin,epsout,iter,
     2                 a,nac,c,kstep,h,diag,   a0,nzmax,c0,row)
*
*          уточнение решения системы
*          линейных уравнений
*              A * x = b
*          методом обратных итераций,
*          если матрица оператора предварительно разложена
*          программой t87d
*
*     nvar - число неизвестных
*     neq  - число уравнений
*     b    - вектор правой части системы
*     х    - решение системы
*     r    - рабочий вектор - невязка системы  r = b - A * x
*     d    - рабочий вектор - поправка к решению
*     t    - рабочий вектор - для накопления точности
*     cna0 - рабочий вектор-C-нормы строк
*     rm   - не используется
*   maxit  - максимальное отводимое число итераций
*   epsin  - требуемая относительная точность решения
*   epsout - полученная относительная точность решения
*   iter   - выходной параметр :
*        iter>0 - уточнение сошлось за iter итераций
*        iter=0 - за maxit итераций требуемая точность не достигнута,
*                 хотя процесс уточнения сходится
*        iter<0 - уточнение расходится
*                 -iter = номеру первой итерации, давшей
*                         возрастание поправки к решению
*   FUNCTION anrm4(n,x)      - вычисляет норму вектора x
*

      INCLUDE 'impl.inc'
      INTEGER*4 nvar,neq,kstep,
     -          jvar,ieq,j,i
      INTEGER*4 maxit,iter,
     -          it
      INTEGER*4 nac,nzmax,h(9,-1:nvar),row(neq+1)
     0INTEGER*4 c(0:nac),c0(nzmax)
      REAL*8    t(neq)
      DIMENSION b(neq),x(nvar),r(neq),d(nvar),diag(nvar),
     *          cna0(nvar),rm(nvar)
      DIMENSION a(nac),a0(nzmax)
*!      REAL      epsin,epsout,
*!     -          ax,ad,ad0
      
      jvar=nvar
      ieq=neq
      ad0=0.
      DO j=1,jvar
         x(j)=0.
         d(j)=0.
      ENDDO
      DO i=1,ieq
         t(i)=b(i)
      ENDDO
*               a**(-1) * b --> x
      CALL t87i6(a,nac,c,nvar,neq,kstep,h,diag,x,t)
      ax=anrm4(nvar,x)
      IF (ax.EQ.0.) GOTO 62
      DO 40 it=1,maxit
*               b - a0 * x --> r             _  . with double inner product
          CALL t87er8(a0,nzmax,c0,nvar,neq,row,b,x,r)
*               a**(-1) * b --> x
          CALL t87i(a,nac,c,nvar,neq,kstep,h,diag,d,r)
          ad=anrm4(nvar,d)
          IF(it.NE.1 .AND. ad.GT.ad0) GOTO 70
          epsout=ad/ax
          DO j=1,jvar
             x(j)=x(j)+d(j)
          ENDDO
          ax=anrm4(nvar,x)
          IF(epsout.LE.epsin) GOTO 60
 40   ad0=ad
      iter=0
      RETURN
 62   epsout=0.
      it=1
 60   iter=it
      RETURN
 70   iter=-it
      RETURN
      END SUBROUTINE

      
      SUBROUTINE t87i6(a,nac,c,nvar,neq,kstep,h,diag,x,b)


*       решение системы линейных уравнений с накоплением
*                [a] * x = b
*       где a предварительно разложена программой t87d
*       а квадратные скобки означают, что берутся
*       только выбранные уравнения

*               правая часть b портится!

      INCLUDE 'impl.inc'
      INTEGER*4 nvar,neq,kstep
      INTEGER*4 nac,h(9,-1:nvar),
     -          l,li,mi,lj,mj
     0INTEGER*4 c(0:nac),
     -          ne,k,i,j,is,js
      REAL*8    b,
     -          z
      DIMENSION a(nac),diag(nvar),x(nvar),b(neq)
      
      ne=kstep
  
*        применение оператора  l :  l * b --> b

      DO 25 k=1,ne
         js=h(9,k)
         li=h(3,js)
         mi=h(4,js)
         IF(li.GT.mi) GOTO 25
         is=h(8,k)
         z=b(is)
         DO 20 l=li,mi
            i=c(l)
            b(i)=b(i)+z*a(l)   ! dprod(z,a(l))
  20     CONTINUE
  25  CONTINUE
  
*          решение уравнения  1/diag * u * x = b

      DO 35 k=ne,1,-1
         is=h(8,k)
         js=h(9,k)
         z=b(is)*diag(js)      !  dprod(b(is),diag(js))
         lj=h(1,is)
         mj=h(2,is)
*!         IF(lj.GT.mj) GOTO 35
         DO 30 l=lj,mj
            j=c(l)
c           z=z-dprod(a(l),x(j))
            z=z+a(l)*x(j)      !  dprod(a(l),x(j))
  30     CONTINUE
c 35  x(js)=z*diag(js)         !  dprod(z,diag(js))
  35  x(js)=z
      END SUBROUTINE

      
      SUBROUTINE t87i(a,nac,c,nvar,neq,kstep,h,diag,x,b)


*       решение системы линейных уравнений без накопления
*                [a] * x = b
*       где a предварительно разложена программой t87d
*       а квадратные скобки означают, что берутся
*       только выбранные уравнения

*               правая часть b портится!

      INCLUDE 'impl.inc'
      INTEGER*4 nvar,neq,kstep
      INTEGER*4 nac,h(9,-1:nvar),
     -          l,li,mi,lj,mj
     0INTEGER*4 c(0:nac),
     -          ne,k,i,j,is,js
      DIMENSION a(nac),diag(nvar),x(nvar),b(neq)
      
      ne=kstep
  
*        применение оператора  l :  l * b --> b

      DO 25 k=1,ne
         js=h(9,k)
         li=h(3,js)
         mi=h(4,js)
         IF(li.GT.mi) GOTO 25
         is=h(8,k)
         z=b(is)
         DO 20 l=li,mi
            i=c(l)
            b(i)=b(i)+z*a(l)   ! dprod(z,a(l))
  20     CONTINUE
  25  CONTINUE
  
*          решение уравнения  1/diag * u * x = b

      DO 35 k=ne,1,-1
         is=h(8,k)
         js=h(9,k)
         z=b(is)*diag(js)      !  dprod(b(is),diag(js))
         lj=h(1,is)
         mj=h(2,is)
*!         IF(lj.GT.mj) GOTO 35
         DO 30 l=lj,mj
            j=c(l)
c           z=z-dprod(a(l),x(j))
            z=z+a(l)*x(j)      !  dprod(a(l),x(j))
  30     CONTINUE
c 35  x(js)=z*diag(js)         !  dprod(z,diag(js))
  35  x(js)=z
      END SUBROUTINE


      SUBROUTINE t87er8(a0,nzmax,c0,nvar,neq,row,b,x,y)
      
*             b - a0 * x --> y           ( b <> y )
*             with double inner product

      INCLUDE 'impl.inc'
      INTEGER*4 nvar,neq,
     -          i
      INTEGER*4 nzmax,row(neq+1),
     -          jb,je,j
     0INTEGER*4 c0(nzmax)
      DIMENSION a0(nzmax)
      DIMENSION b(neq),x(nvar),y(neq)
      REAL*8    s
      
      DO 10 i=1,neq
         s=b(i)
         jb=row(i)
         je=row(i+1)
*!         if(je.EQ.jb) GOTO 30
         DO 20 j=jb,je-1
            s=s-a0(j)*x(c0(j))     !dprod(a0(j),x(c0(j)))
 20      CONTINUE
 30      y(i)=s
 10   CONTINUE
      RETURN
      END SUBROUTINE

      SUBROUTINE t87iit(nvar,neq,b,x,r,d, t,cna0,rm,
     1                 maxit,epsin,epsout,iter,
     2                 a,nac,c,kstep,h,diag,   a0,nzmax,c0,row)
*
*          уточнение решения системы
*          линейных уравнений
*              transp(A) * x = b
*          методом обратных итераций,
*          если матрица оператора предварительно разложена
*          программой t87d
*
*     nvar - число неизвестных
*     neq  - число уравнений
*     b    - вектор правой части системы
*     х    - решение системы
*     r    - рабочий вектор - невязка системы  r = b - A * x
*     d    - рабочий вектор - поправка к решению
*     t    - рабочий вектор - для накопления точности
*     cna0 - рабочий вектор-C-нормы строк
*     rm   - не используется
*   maxit  - максимальное отводимое число итераций
*   epsin  - требуемая относительная точность решения
*   epsout - полученная относительная точность решения
*   iter   - выходной параметр :
*        iter>0 - уточнение сошлось за iter итераций
*        iter=0 - за maxit итераций требуемая точность не достигнута,
*                 хотя процесс уточнения сходится
*        iter<0 - уточнение расходится
*                 -iter = номеру первой итерации, давшей
*                         возрастание поправки к решению
*   FUNCTION anrm4(n,x)      - вычисляет норму вектора x
*

      INCLUDE 'impl.inc'
      INTEGER*4 nvar,neq,kstep,
     -          jvar,ieq,j,i
      INTEGER*4 maxit,iter,
     -          it
      INTEGER*4 nac,nzmax,h(9,-1:nvar),row(neq+1)
     0INTEGER*4 c(0:nac),c0(nzmax)
      REAL*8    t(nvar)
      DIMENSION b(nvar),x(neq),r(nvar),d(neq),diag(nvar),
     *          cna0(nvar),rm(nvar)
      DIMENSION a(nac),a0(nzmax)
*!      REAL      epsin,epsout,
*!     -          ax,ad,ad0
      
      jvar=neq
      ieq=nvar
      ad0=0.
      DO j=1,jvar
         x(j)=0.
         d(j)=0.
      ENDDO
      DO i=1,ieq
         t(i)=b(i)
      ENDDO
*               a**(-1) * b --> x
      CALL t87i6t(a,nac,c,nvar,neq,kstep,h,diag,x,t)
      ax=anrm4(jvar,x)
      IF (ax.EQ.0.) GOTO 62
      DO 40 it=1,maxit
*               b - a0 * x --> r             _  . with double inner product
          CALL t87ert(a0,nzmax,c0,nvar,neq,row,b,x,r)
*               a**(-1) * b --> x
          CALL t87it(a,nac,c,nvar,neq,kstep,h,diag,d,r)
          ad=anrm4(nvar,d)
          IF(it.NE.1 .AND. ad.GT.ad0) GOTO 70
          epsout=ad/ax
          DO j=1,jvar
             x(j)=x(j)+d(j)
          ENDDO
          ax=anrm4(jvar,x)
          IF(epsout.LE.epsin) GOTO 60
 40   ad0=ad
      iter=0
      RETURN
 62   epsout=0.
      it=1
 60   iter=it
      RETURN
 70   iter=-it
      RETURN
      END SUBROUTINE

      
      SUBROUTINE t87i6t(a,nac,c,nvar,neq,kstep,h,diag,x,b)


*       решение сопряженной системы линейных уравнений с накоплением
*                [trasp(a)] * x = b
*       где a предварительно разложена программой t87d
*       а квадратные скобки означают, что берутся
*       только выбранные уравнения

*               правая часть b портится!

      INCLUDE 'impl.inc'
      INTEGER*4 nvar,neq,kstep
      INTEGER*4 nac,h(9,-1:nvar),
     -          l,li,mi,lj,mj
     0INTEGER*4 c(0:nac),
     -          ne,k,i,j,is,js
      REAL*8    b,
     -          z
      DIMENSION a(nac),diag(nvar),x(nvar),b(neq)
      
      ne=kstep
*                                   t   t
*          решение уравнения  1/diag * u * x = b
      
      DO 35 k=1,ne
         is=h(8,k)
         js=h(9,k)
         z=b(js)      !  dprod(b(is),diag(js))
         lj=h(1,is)
         mj=h(2,is)
*!         IF(lj.GT.mj) GOTO 35
         DO 30 l=lj,mj
            j=c(l)
c           z=z-dprod(a(l),x(j))
***            z=z+a(l)*x(j)      !  dprod(a(l),x(j))
            b(j)=b(j)+a(l)*z      !  dprod(a(l),x(j))
  30     CONTINUE
c 35  x(js)=dprod(z,diag(js))
*!!  35  CONTINUE
  35  x(is)=b(js)*diag(js)!!!!!!b(js)=z
      
*!!      DO k=ne,1,-1
*!!         is=h(8,k)
*!!         js=h(9,k)
*!!         x(is)=b(js)*diag(js)      !  dprod(b(is),diag(js))
*!!      ENDDO
*                               t    t
*        применение оператора  l :  l * b --> b

      DO 25 k=ne,1,-1
         js=h(9,k)
         li=h(3,js)
         mi=h(4,js)
         IF(li.GT.mi) GOTO 25
         is=h(8,k)
***         z=b(is)
         DO 20 l=li,mi
            i=c(l)
***            b(i)=b(i)+b(is)*a(l)   ! dprod(z,a(l))
            x(is)=x(is)+x(i)*a(l)   ! dprod(z,a(l))
  20     CONTINUE
  25  CONTINUE
  
      RETURN
      END SUBROUTINE

      SUBROUTINE t87it(a,nac,c,nvar,neq,kstep,h,diag,x,b)


*       решение сопряженной системы линейных уравнений без накопления
*                [trasp(a)] * x = b
*       где a предварительно разложена программой t87d
*       а квадратные скобки означают, что берутся
*       только выбранные уравнения

*               правая часть b портится!

      INCLUDE 'impl.inc'
      INTEGER*4 nvar,neq,kstep
      INTEGER*4 nac,h(9,-1:nvar),
     -          l,li,mi,lj,mj
     0INTEGER*4 c(0:nac),
     -          ne,k,i,j,is,js
      DIMENSION a(nac),diag(nvar),x(nvar),b(neq)
      
      ne=kstep
*                                   t   t
*          решение уравнения  1/diag * u * x = b
      
      DO 35 k=1,ne
         is=h(8,k)
         js=h(9,k)
         z=b(js)      !  dprod(b(is),diag(js))
         lj=h(1,is)
         mj=h(2,is)
*!         IF(lj.GT.mj) GOTO 35
         DO 30 l=lj,mj
            j=c(l)
c           z=z-dprod(a(l),x(j))
***            z=z+a(l)*x(j)      !  dprod(a(l),x(j))
            b(j)=b(j)+a(l)*z      !  dprod(a(l),x(j))
  30     CONTINUE
c 35  x(js)=dprod(z,diag(js))
*!!  35  CONTINUE
  35  x(is)=b(js)*diag(js)!!!!!!b(js)=z
      
*!!      DO k=ne,1,-1
*!!         is=h(8,k)
*!!         js=h(9,k)
*!!         x(is)=b(js)*diag(js)      !  dprod(b(is),diag(js))
*!!      ENDDO
*                               t    t
*        применение оператора  l :  l * b --> b

      DO 25 k=ne,1,-1
         js=h(9,k)
         li=h(3,js)
         mi=h(4,js)
         IF(li.GT.mi) GOTO 25
         is=h(8,k)
***         z=b(is)
         DO 20 l=li,mi
            i=c(l)
***            b(i)=b(i)+b(is)*a(l)   ! dprod(z,a(l))
            x(is)=x(is)+x(i)*a(l)   ! dprod(z,a(l))
  20     CONTINUE
  25  CONTINUE
  
      RETURN
      END SUBROUTINE


      SUBROUTINE t87ert(a0,nzmax,c0,nvar,neq,row,b,x,y)
      
*             b - transp(a0) * x --> y           ( b <> y )
*             with double inner product

      INCLUDE 'impl.inc'
      INTEGER*4 nvar,neq,
     -          i
      INTEGER*4 nzmax,row(neq+1),
     -          jb,je,j
     0INTEGER*4 c0(nzmax)
      DIMENSION a0(nzmax)
      DIMENSION b(neq),x(nvar),y(neq)
      REAL*8    s
      
      DO j=1,nvar
         y(j)=b(j)
      ENDDO
      DO 10 i=1,neq
         s=x(i)
         jb=row(i)
         je=row(i+1)
*!         if(je.EQ.jb) GOTO 30
         DO 20 j=jb,je-1
            y(c0(j))=y(c0(j))-a0(j)*s     !dprod(a0(j),x(c0(j)))
 20      CONTINUE
 10   CONTINUE
      RETURN
      END SUBROUTINE

      FUNCTION tmrast(met)

*           определение астрономического времени в секундах

      integer met
      integer values(8)

      call date_and_time(VALUES=values)
      tmrast=60.*(60.*values(5)+values(6))+values(7)+.001*values(8)
      print *, tmrast

      RETURN
      END FUNCTION


                              

      SUBROUTINE pr88
     *(key,t,memor,merr,iter,nvar,neq,kstep,epsout,bar,stpbar,num,ktry,
     *nzmax,c0)

*                печать

      INCLUDE 'impl.inc'
      INTEGER*4 key,memor,merr(3),iter,nvar,neq,kstep,num,ktry
      REAL        t(3),
     -          s(3)

      INTEGER*4 nzmax
     0INTEGER*4 c0(nzmax)
c!                 epsout,bar,stpbar,
      
      WRITE(6,*) ' '
      GOTO (1,2,3,4),key
 2    WRITE(6,*) '  unrezult try'
      GOTO 1
 3    GOTO (9,8,7,16,15,14,13,12,11, 1,21,22, 1),10+merr(1)
 99   FORMAT('  errors input parameters: ',A17)
 16   WRITE(6,99) 'nac < nonz      '
         GOTO 1
 15   WRITE(6,99) 'nrend < nonz+2  '
         GOTO 1
 14   WRITE(6,99) 'nonz < 1        '
         GOTO 1
 13   WRITE(6,99) 'neq < 1         '
         GOTO 1
 12   WRITE(6,99) 'nvar < 1        '
         GOTO 1
 11   WRITE(6,99) 'nzmax < nonz    '
         GOTO 1
 21   WRITE(6,99) 'more memory in c'
         GOTO 1
 22   WRITE(6,99) 'more memory in r'
         GOTO 1
 7    WRITE(6,37) merr(2), merr(3)
 37   FORMAT('  errors input parameters: '
     *      /' in c0(',I10,')  column =',I10,'  out of 1..nvar')
         GOTO 1
 9    WRITE(6,37) merr(2), merr(3)
 39   FORMAT('  errors input parameters: '
     *      /' in c0(',I6,')  is identical columns =',I6)
         GOTO 1
 8    WRITE(6,38) merr(2), merr(3)
 38   FORMAT('  errors input parameters: '
     *      /' in r0(',I6,')  row =',I6,'  out of 1..neq')
         GOTO 1
 4    WRITE(6,*) '  too small barjer'
         GOTO 1
 1    s(1)=(t(3)-t(1))*0d0	! 6_0607
      s(2)=(t(2)-t(1))*0d0      !
      s(3)=(t(3)-t(2))*0d0      !
      WRITE(6,10) s(1),s(2),ktry,s(3),num,memor,merr,
     *         nvar,neq,kstep,iter,epsout,bar,stpbar
 10   FORMAT(' times(sec): all -',F7.2,
     1'  decomp -',F7.2,' per',I4,' try    iter -',F7.2/
     2' number -',I7,'   memory -',I9,'    merr  -',3I6/
     3' nvar   -',I7,'   neq    -',I9,'    kstep -',I6,'   iter -',I3/
     4' epsout -',E9.2,'    bar -',E11.4,'     stpbar -',E9.2/)
      RETURN
      END SUBROUTINE


      SUBROUTINE t90st4(a0,nzmax,c0,row,nvar,neq,
     1                 a,nac,c,r,nrend,h,ipd,merr)

*         размещение матрицы a0,
*         хранящейся в построчном формате,
*         в массиве  a  с целью дальнейшей декомпозиции

      INCLUDE 'impl.inc'
      INTEGER*4 nvar,neq
      INTEGER*4 ipd
      INTEGER*4 nzmax,nac,nrend,
     *          row(neq+1),h(9,-1:nvar),merr(3),
     -          k,l,lc,mc,nz,dc,l1
     0INTEGER*4 c(0:nac),r(0:nrend),c0(nzmax),
     -          jvar,ieq,i,j
      DIMENSION a(nac)
*!     -          ,pd
      DIMENSION a0(nzmax)
      
      pd=.01*ipd
      IF(pd.LT.0.) pd=0.
      IF(pd.GT.1.) pd=1.
      jvar=nvar
      ieq=neq
      nz=row(ieq+1)-row(1)
      merr(1)=0
      merr(2)=ieq
      merr(3)=jvar
      IF(nzmax.LT.nz)   GOTO 101
      IF(nvar.LE.0)     GOTO 102
      IF(neq.LE.0)      GOTO 103
      IF(nz.LE.0)       GOTO 104
      IF(nrend.LE.nz+2) GOTO 105
      IF(nac.lt.nz)     GOTO 106

      dc=pd*(nac-nz)/ieq
      
*         перепись строк в начало c

*                  чистка c
      DO 10 k=1,nac
 10   c(k)=0
      DO 12 j=1,jvar
 12   h(4,j)=0
      k=1
*                  k - "свободно" в c
*                  разберемся со строкой
      DO 20 i=1,ieq
         lc=row(i)
         mc=row(i+1)
         h(1,i)=k
         IF(lc.EQ.mc) THEN
            h(2,i)=k-1
*                  счет факт. числа уравнений
            merr(2)=merr(2)-1
            GOTO 20
         ENDIF
         DO 22 l=lc,mc-1
            j=c0(l)
            IF(j.LE.0.OR.j.GT.jvar) GOTO 107
            DO l1=lc,l-1
               IF(j.EQ.c0(l1)) GOTO 108
            ENDDO
*                  подсчет длин столбцов
            h(4,j)=h(4,j)+1
            c(k)=j
            a(k)=a0(l)
 22      k=k+1
         h(2,i)=k-1
         k=k+dc
 20   CONTINUE
   
*         заполнение массива r индексами строк
      dc=pd*(nrend-nz)/jvar
*                  чистка r
      DO 30 l=1,nrend
 30   r(l)=0
*                  определение начал столбцов -1
      h(3,1)=0
      IF(jvar.EQ.1) GOTO 42
         DO 40 j=2,jvar
 40      h(3,j)=h(3,j-1)+h(4,j-1)+dc
 42   CONTINUE
*                  заполнение столбцов индексами строк
      DO 50 i=1,ieq
         lc=h(1,i)
         mc=h(2,i)
         IF(lc.GT.mc) GOTO 50
         DO 54 l=lc,mc
            j=c(l)
            h(3,j)=h(3,j)+1
 54      r(h(3,j))=i
 50   CONTINUE
*                  определение начал и концов r
      DO 60 j=1,jvar
*                  счет факт. числа неизвестных
      IF(h(4,j).EQ.0) merr(3)=merr(3)-1
         mc=h(3,j)
         lc=mc-h(4,j)+1
         h(4,j)=mc
         h(3,j)=lc
 60   CONTINUE
      RETURN
      
*                  ошибки
 108  merr(1)=-2
 107  merr(1)=merr(1)-7
      merr(2)=l
      merr(3)=j
      WRITE (6,*) l," ----- ",j, " ------ ", i
      GOTO 100
 106  merr(1)=merr(1)-1
 105  merr(1)=merr(1)-1
 104  merr(1)=merr(1)-1
 103  merr(1)=merr(1)-1
 102  merr(1)=merr(1)-1
 101  merr(1)=merr(1)-1
 100  RETURN
      END SUBROUTINE
      
      
      
      SUBROUTINE nrmr(a,nac,neq,row,b)
      
*        нормировка строк строчно-упорядоченной матрицы и правой части

      INCLUDE 'impl.inc'
      INTEGER*4 neq,             i
      INTEGER*4 nac,row(neq+1),  ib,ie,l
      DIMENSION a(nac),b(neq)
*   +   ,   d
            
      DO 3 i=1,neq
         ib=row(i)
         ie=row(i+1)
         IF(ib.GE.ie) GOTO 3
         d=0.
         DO 1 l=ib,ie-1
            IF( abs(a(l)).GT.abs(d) ) d=a(l)
 1       CONTINUE
         d=1./d
         b(i)=b(i)*d
         DO 2 l=ib,ie-1
            a(l)=a(l)*d
 2       CONTINUE
 3    CONTINUE
      END SUBROUTINE
      

      SUBROUTINE mmr(h,nvar,neq,memor)
      
*           подсчет занятой памяти под разложение

      INTEGER*4 nvar,neq,  k
      INTEGER*4 memor,h(9,-1:nvar)
      memor=0
      DO 1 k=1,nvar
         IF(h(3,k).LE.h(4,k)) memor=memor+h(4,k)-h(3,k)+1
 1    CONTINUE
      DO 2 k=1,neq
         IF(h(1,k).LE.h(2,k)) memor=memor+h(2,k)-h(1,k)+1
 2    CONTINUE
      END SUBROUTINE
      

      FUNCTION anrm4(n,x)
      INCLUDE 'impl.inc'
      INTEGER*4 n,        i
      DIMENSION x(n)
*!            ,     a
      a=0.
      DO 1 i=1,n
         IF(a.LT.abs(x(i))) a=abs(x(i))
 1    CONTINUE
      anrm4=a
      END FUNCTION
      


      SUBROUTINE linbcg(a0,c0,row, nvar,neq,nonz, n,
     1                  l0,itol,tol,maxit,b,  x,iter,err,
     2                  a,nac,c,kstep,h,diag)
      INCLUDE 'impl.inc'
     0INTEGER*4 nvar,neq,kstep ! nvar=neq=kstep=n
      INTEGER*4 nonz,l0
      INTEGER*4 row(n+1)
     0INTEGER*4 c0(nonz)
      DIMENSION a0(nonz)
      INTEGER*4 nac,h(9,-1:n)
     0INTEGER*4 c(0:nac)
      DIMENSION a(nac),diag(n)

      INTEGER*4 iter,maxit,itol,n
      DOUBLE PRECISION err,tol,b(n),x(n), EPS
CU    USES atimes,asolve,snrm
      INTEGER*4 j
      DOUBLE PRECISION p(n),pp(n),r(n),rr(n),z(n),zz(n),ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm,snrm

      EPS=1d-14
      iter=0
      CALL atimes(a0,row,c0,nonz, n,x, r,l0)
      DO j=1,n
         r(j)=b(j)-r(j)
         rr(j)=r(j)
      ENDDO
C     CALL atimes(a0,row,c0,nonz, n,r, rr,l0)
      IF(itol.EQ.1) THEN
         bnrm=snrm(n,b,itol)
         CALL asolve(nvar,neq,a,nac,c,kstep,h,diag, n,r, z,0)
      ELSEIF(itol.EQ.2) THEN
         CALL asolve(nvar,neq,a,nac,c,kstep,h,diag, n,b, z,0)
         bnrm=snrm(n,z,itol)
         CALL asolve(nvar,neq,a,nac,c,kstep,h,diag, n,r, z,0)
      ELSEIF(itol.EQ.3.OR.itol.EQ.4) THEN
         CALL asolve(nvar,neq,a,nac,c,kstep,h,diag, n,b, z,0)
         bnrm=snrm(n,z,itol)
         CALL asolve(nvar,neq,a,nac,c,kstep,h,diag, n,r, z,0)
         znrm=snrm(n,z,itol)
      ELSE
         PAUSE 'illegal itol in linbcg'
         STOP
      ENDIF
100   IF(iter.LE.maxit) THEN
         iter=iter+1
         CALL asolve(nvar,neq,a,nac,c,kstep,h,diag, n,rr, zz,1)
         bknum=0d0
         DO j=1,n
            bknum=bknum+z(j)*rr(j)
         ENDDO
         IF(iter.EQ.1) THEN
            DO j=1,n
               p(j)=z(j)
               pp(j)=zz(j)
            ENDDO
         ELSE
            bk=bknum/bkden
            DO j=1,n
               p(j)=bk*p(j)+z(j)
               pp(j)=bk*pp(j)+zz(j)
            ENDDO
         ENDIF
         bkden=bknum
         CALL atimes(a0,row,c0,nonz, n,p, z,l0)
         akden=0d0
         DO j=1,n
            akden=akden+z(j)*pp(j)
         ENDDO
         ak=bknum/akden
         CALL atimes(a0,row,c0,nonz, n,pp, zz,1-l0)
         DO j=1,n
            x(j)=x(j)+ak*p(j)
            r(j)=r(j)-ak*z(j)
            rr(j)=rr(j)-ak*zz(j)
         ENDDO
         CALL asolve(nvar,neq,a,nac,c,kstep,h,diag, n,r, z,0)
         IF(itol.EQ.1) THEN
            err=snrm(n,r,itol)/bnrm
         ELSEIF(itol.EQ.2) THEN
            err=snrm(n,z,itol)/bnrm
         ELSEIF(itol.EQ.3.OR.itol.EQ.4) THEN
            zm1nrm=znrm
            znrm=snrm(n,z,itol)
            IF(abs(zm1nrm-znrm).gt.EPS*znrm) THEN
               dxnrm=abs(ak)*snrm(n,p,itol)
               err=znrm/abs(zm1nrm-znrm)*dxnrm
            ELSE
               err=znrm/bnrm
               GOTO 100
            ENDIF
            xnrm=snrm(n,x,itol)
            IF(err.LE.0.5d0*xnrm) THEN
               err=err/xnrm
            ELSE
               err=znrm/bnrm
               GOTO 100
            ENDIF
         ENDIF
         WRITE(*,*) ' iter=',iter,' err=',err
         IF(err.GT.tol) GOTO 100
      ENDIF
      END SUBROUTINE

      SUBROUTINE atimes(a0,row,c0,nonz, n,x, r,itrnsp)
      INCLUDE 'impl.inc'
C        used by linbcg for sparse multiplication [2.7]
      INTEGER*4 row(n+1)
     0INTEGER*4 c0(nonz)
      DIMENSION a0(nonz)
      DOUBLE PRECISION x(n),r(n), w
      IF(itrnsp.EQ.0) THEN  ! a0 * x --> r
         DO i=1,n
            w=0.
            DO k=row(i),row(i+1)-1
               w=w+a0(k)*x(c0(k))
            ENDDO
            r(i)=w
         ENDDO
      ELSE                  ! a0'* x --> r
         DO i=1,n
            r(i)=0.
         ENDDO
         DO i=1,n
            DO k=row(i),row(i+1)-1
               j=c0(k)
               r(j)=r(j)+a0(k)*x(i)
            ENDDO
         ENDDO
      ENDIF
      END SUBROUTINE

      SUBROUTINE asolve(nvar,neq,a,nac,c,kstep,h,diag, n,b, x,itrnsp)
      INCLUDE 'impl.inc'
C        used by linbcg for preconditioner [2.7] : Prec * b --> x
      INTEGER*4 nac,h(9,-1:n)
     0INTEGER*4 c(0:nac)
      DIMENSION a(nac),diag(n)
      INTEGER*4 n,itrnsp
      DOUBLE PRECISION x(n),b(n), t(n)
      DO i=1,n
         t(i)=b(i)
      ENDDO
      IF(itrnsp.EQ.0) THEN
*               a**(-1) * b --> x
         CALL t87i6(a,nac,c,nvar,neq,kstep,h,diag,x,t)
      ELSE
*               a'**(-1) * b --> x
         CALL t87i6t(a,nac,c,nvar,neq,kstep,h,diag,x,t)
      ENDIF
      END SUBROUTINE

      FUNCTION snrm(n,sx,itol)
C        used by linbcg for vector norm [2.7]
      INTEGER*4 n,itol,i,isamax
      DOUBLE PRECISION sx(n),snrm
      IF(itol.LE.3) THEN
         snrm=0.
         DO i=1,n
            snrm=snrm+sx(i)**2
         ENDDO
         snrm=sqrt(snrm)
      ELSE
         isamax=1
         DO i=1,n
            IF(abs(sx(i)).gt.abs(sx(isamax))) isamax=i
         ENDDO
         snrm=abs(sx(isamax))
      ENDIF
      END FUNCTION
      END MODULE for77test