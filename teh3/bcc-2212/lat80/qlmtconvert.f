C--------------------------------------- 
C This program converts the LAPW qlmt file, C to run the dosapwn.f program
C------------------------------------------
        program convert
        implicit real*8 (a-h,o-z)
        open(1,file='qlmt')
        open(2,file='dosapw.eng')
        read(1,*)  NB
        read(1,*) C 
 100    continue
        read(1,*,end=200) x,x,x,x,nband
        NB=9
        if (nband .lt. NB) then
        PRINT*, 'YOU HAVE LESS THAN ',NB,' BANDS IN YOUR
     C FILE'
        endif
        do nn=1,nband
        read(1,*) eign
        read(1,*) s1,p1,d1,d11
        read(1,*) s2,p2,d2,d22
c       read(1,*) s3,p3,d3,d33
c       read(1,*) s4,p4,d4,d44
c       read(1,*) s5,p5,d5,d55
        if (nn .le. NB) then
           write(2,'(6F8.5)')eign,s1,p1,d1,d11,0.0
           write(2,'(6F8.5)')eign,s2,p2,d2,d22,0.0
c          write(2,'(6F8.5)')eign,s3,p3,d3,d33,0.0
c          write(2,'(6F8.5)')eign,s4,p4,d4,d44,0.0
c          write(2,'(6F8.5)')eign,s5,p5,d5,d55,0.0
        endif
        enddo
        go to 100
  200 continue
        end

