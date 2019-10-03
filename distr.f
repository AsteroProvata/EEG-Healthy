         program distribution
c---- reads the distributed quantity  xm(i).  
c--- Finds the size distribution of  xm(i).
                 
        implicit real*8 (a-h, o-z)
        integer*4 seed

        parameter(L=80,L2=L*L)
        parameter (idistr=500,distr=idistr)
        dimension xm(L2),xseg(idistr)
        dimension px(idistr)
         
c        dimension b(maxconf,maxconf),ssk(maxconf)
c        dimension wclustup(maxconf),wclustdown(maxconf),wclust(maxconf)
c        dimension lfirst(kwindow),lsecond(kwindow)
c        dimension bp(0:idistr),psk(idistr),pwc(idistr)


        open(unit=11, file='00438_40Hz_matrix.txt', status='old')
        open(unit=99, file='./Distributions/distr00438_40Hz_matri0.dat')
       
c***************Initiation des variables ***************
c              do i=1,2
c              read(11,*)
c              enddo
              do i=1,L2
              xm(i)=0.0
              enddo
c
              do i=1,idistr
              px(i)=0.0
              enddo



c-------------Read the data file---------------------
            
                ra=rand(seed)
                xtot=0.0
                do  k=1,L2
                read(11,*,end=18) i,j,xm(k) 
                xtot=xtot+xm(k)
cc                 write(*,*) i,j,xm(k)
                enddo
              
c**for random***
c                do i=1,L
c                 do j=1,L
c                ra=rand()
c                re=rand()
c                rm(i,j)=ra
c                xm(i,j)=re
c                enddo
c                enddo
18              Lmax=k-1
           



c--------------do distributions of xm(i)-------------------------------------

               xmax=-10000000
               xmin=10000000
               do i=1,L2
               if(xm(i).gt.xmax) xmax=xm(i)
               if(xm(i).lt.xmin) xmin=xm(i)
               enddo
               write(*,*) 'help1'
c
               xbin=(xmax-xmin)/distr
c
               xseg(1)=xmin
               do i=2,idistr
               xseg(i)=xseg(i-1)+xbin
               enddo
               write(*,*) 'help2',xmin,xmax
c
               do i=1,L2
               do k=1,idistr-1
                  if(xm(i).ge.xseg(k).and.xm(i).lt.xseg(k+1)) then
                       px(k)=px(k)+1.
                       go to 38
                    endif
38             enddo
               enddo
               write(*,*) 'help3'
               

               xtotal=0.0
               do k=1,idistr
               xtotal=xtotal+xseg(k)
               enddo 
               write(*,*) 'xtotal=',xtotal,xtot

               do k=1,idistr
               px(k)=px(k)/dble(Lmax)
               enddo 
               write(*,*) 'matrix size distribution done'
               do k=2,idistr
               if(px(k).ne.0) write(99,*) xseg(k),px(k),k
               enddo

                            
c-----------calculate cummulative-----------------------------------------------------------
c               cum=0.0
c               do k=idistr,1,-1
c                cum=cum+px(k)
c                if(cum.gt.0.0) write(99,*) xseg(k),px(k),cum,k
c               enddo
                           
c----------------------------------------------------------------------------
c---------- Calculate degree(weight of links) on each vertex-----------------
c
c
c
c               do i=1,L
c               do j=1,L
c               dr(i)=dr(i)+rm(i,j)
c               dx(i)=dx(i)+xm(i,j)
c               enddo
c               write(77,*) i,dr(i),dx(i)
c               enddo



           stop
           end


