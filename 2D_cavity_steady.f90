!��������������������������������������xxx
!������������㷨����ӭ������ɢ��ʽ��
!�˹�ѹ���㷨��ⷽǻ�������⣨���԰汾��
!����������������������������������
      program QUICK_cavity
      parameter(mx=101,my=101)
      implicit double precision(a-h,o-z)
	  double precision u(mx,my+1),v(mx+1,my),p(mx+1,my+1)
	  double precision  un(mx,my+1),vn(mx+1,my),pn(mx+1,my+1)
      common /ini/u,v,p
      c2=2.25 
      re=100.0
      dt=0.0005
      dx=1.0/float(mx-1)
      dy=1.0/float(my-1)
!----------------------------------------------------------------------------------------
!  u��v��pΪtʱ��ֵ��un��vn��pnΪt+1ʱ��ֵ��
!  mx��myΪ�����������c2Ϊ����ѹ��ϵ����ƽ����reΪ��ŵ����
!-----------------------------------------------------------------------------------------
      num=0
      err=100.00
!nun����������err���ж��˹�ѹ������������ı�׼
      call initial
!�����ʼ����������Ϊ�˹�ѹ���㷨���
      do while(err.gt.1e-4.and.num.lt.1e6)
         err=0.0
	     call quick(u,v,p,un,vn,mx,my,dx,dy,dt,re)
!QUICK��ɢ��ʽ��⶯�����̣��õ�un��vn
	     call calp(p,un,vn,pn,mx,my,dx,dy,dt,c2)
!��ѹǿpn
         call check(u,v,p,un,vn,pn,mx,my,dx,dy,dt,c2,err)
!У��������Ϣ���ж��Ƿ�������ͬʱ����u��v��p
	     write(*,*) 'error=',err
	     num=num+1
         write(*,*) num
!��Ļ�������
      enddo
      call output(u,v,p,mx,my,dx,dy)
!�������ļ�
      End
!
!
      subroutine initial
!��ʼ������
	  parameter(mx=101,my=101)
	  double precision u(mx,my+1),v(mx+1,my),p(mx+1,my+1)
	  common /ini/u,v,p
      do i=1,mx+1
	     do j=1,my+1
	        p(i,j)=1.0
	     enddo
	  enddo
	  do i=1,mx
	     do j=1,my+1
	        u(i,j)=0.0
		    if(j.eq.my+1) u(i,j)=4.0/3.0
		    if(j.eq.my) u(i,j)=2.0/3.0
	     enddo
	  enddo
	  do i=1,mx+1
	     do j=1,my
		    v(i,j)=0.0
	     enddo
	  enddo
	  endsubroutine
!
!
      subroutine quick(u,v,p,un,vn,mx,my,dx,dy,dt,re)
!��QUICK��ʽ��ɢ��������
	  implicit double precision(a-h,o-z)
	  dimension u(mx,my+1),v(mx+1,my),p(mx+1,my+1),un(mx,my+1),vn(mx+1,my)
      double precision miu
      miu=1.0/re
!�������x�����ٶ�un----------------------------------------------------------------------------
     do i=3,mx-2
         do j=3,my-1
            fw=0.5*(u(i-1,j)+u(i,j))*dy
		    fe=0.5*(u(i,j)+u(i+1,j))*dy
		    fs=0.5*(v(i,j-1)+v(i+1,j-1))*dx
		    fn=0.5*(v(i,j)+v(i+1,j))*dx
		    df=fe-fw+fn-fs
            aw=miu+0.750*alfa(fw)*fw+0.125*alfa(fe)*fe+0.375*(1.0-alfa(fw))*fw
		    aww=-0.125*alfa(fw)*fw
		    ae=miu-0.375*alfa(fe)*fe-0.750*(1.0-alfa(fe))*fe-0.125*(1.0-alfa(fw))*fw
		    aee=0.125*(1.0-alfa(fe))*fe
            as=miu+0.750*alfa(fs)*fs+0.125*alfa(fn)*fn+0.375*(1.0-alfa(fs))*fs
		    ass=-0.125*alfa(fs)*fs
		    an=miu-0.375*alfa(fn)*fn-0.750*(1.0-alfa(fn))*fn-0.125*(1.0-alfa(fs))*fs
		    ann=0.125*(1.0-alfa(fn))*fn
		    ap=aw+ae+as+an+aww+aee+ass+ann+df
!aw��ae��as��an...��Ϊ��������㷨�и���ϵ�������ǰ������ӭ����QUICK��ɢ��ʽ
		    un(i,j)=u(i,j)+dt/dx/dy*(-ap*u(i,j)+aw*u(i-1,j)+ae*u(i+1,j)   &
+as*u(i,j-1)+an*u(i,j+1)+aww*u(i-2,j)+aee*u(i+2,j)  &
+ass*u(i,j-2)+ann*u(i,j+2))-dt*(p(i+1,j)-p(i,j))/dx
	     enddo
      enddo
!-------------------------------------------------------------------------------------------
      j=2
      do i=3,mx-2
         call upbound_u(u,v,p,un,mx,my,dx,dy,dt,re,i,j)
      enddo
	  j=my
      do i=3,mx-2
	     call upbound_u(u,v,p,un,mx,my,dx,dy,dt,re,i,j)
      enddo
      i=2
      do j=2,my
	     call upbound_u(u,v,p,un,mx,my,dx,dy,dt,re,i,j)
      enddo
	  i=mx-1
      do j=2,my
	     call upbound_u(u,v,p,un,mx,my,dx,dy,dt,re,i,j)
      enddo
!�ڲ�߽���һ��ӭ������ɢ��ʽ�õ�----------------------------------------------------
!-------------------------------------------------------------------------------------------
      do i=2,mx-1
         un(i,1)=-un(i,2)
	     un(i,my+1)=2.0-un(i,my)
      enddo
      do j=1,my+1
         un(1,j)=0.0
	     un(mx,j)=0.0
      enddo
!���߽�����������߽���������
!-----------------------------------------------------------------------------
!�������y�����ٶ�vn-----------------------------------------------
      do i=3,mx-1
         do j=3,my-2
	        fw=0.5*(u(i-1,j)+u(i-1,j+1))*dy
		    fe=0.5*(u(i,j)+u(i,j+1))*dy
		    fs=0.5*(v(i,j-1)+v(i,j))*dx
		    fn=0.5*(v(i,j)+v(i,j+1))*dx
		    df=fe-fw+fn-fs
            aw=miu+0.750*alfa(fw)*fw+0.125*alfa(fe)*fe+0.375*(1.0-alfa(fw))*fw
		    aww=-0.125*alfa(fw)*fw
		    ae=miu-0.375*alfa(fe)*fe-0.750*(1.0-alfa(fe))*fe-0.125*(1.0-alfa(fw))*fw
		    aee=0.125*(1.0-alfa(fe))*fe
            as=miu+0.750*alfa(fs)*fs+0.125*alfa(fn)*fn+0.375*(1.0-alfa(fs))*fs
		    ass=-0.125*alfa(fs)*fs
		    an=miu-0.375*alfa(fn)*fn-0.750*(1.0-alfa(fn))*fn-0.125*(1.0-alfa(fs))*fs
		    ann=0.125*(1.0-alfa(fn))*fn
		    ap=aw+ae+as+an+aww+aee+ass+ann+df
		    vn(i,j)=v(i,j)+dt/dx/dy*(-ap*v(i,j)+aw*v(i-1,j)+ae*v(i+1,j)   &
+as*v(i,j-1)+an*v(i,j+1)+aww*v(i-2,j)+aee*v(i+2,j)  &
+ass*v(i,j-2)+ann*v(i,j+2))-dt*(p(i,j+1)-p(i,j))/dy
         enddo
      enddo
!-----------------------------------------------------------------------------
      j=2
	  do i=3,mx-1
         call upbound_v(u,v,p,vn,mx,my,dx,dy,dt,re,i,j)
	  enddo
   	  j=my-1
	  do i=3,mx-1
         call upbound_v(u,v,p,vn,mx,my,dx,dy,dt,re,i,j)
      enddo
      i=2
	  do j=2,my-1
         call upbound_v(u,v,p,vn,mx,my,dx,dy,dt,re,i,j)
	  enddo
	  i=mx
	  do j=2,my-1
         call upbound_v(u,v,p,vn,mx,my,dx,dy,dt,re,i,j)
      enddo
!----------------------------------------------------------------------------
   do i=2,mx
      vn(i,1)=0.0
	  vn(i,my)=0.0
   enddo
   do j=1,my
      vn(1,j)=-vn(2,j)
	  vn(mx+1,j)=-vn(mx,j)
   enddo
!----------------------------------------------------------------------------
   Endsubroutine
!
!
      function alfa(x)
!����,��1��0
      double precision alfa, x
	  if(x.gt.0.d0) then
	    alfa=1.0
	  else
	    alfa=0.0
	  endif
	  end
!
!
      subroutine upbound_u(u,v,p,un,mx,my,dx,dy,dt,re,i,j)
!��һ��ӭ������ɢ��ʽ�õ��ڲ�߽�un��ֵ
	  implicit double precision(a-h,o-z)
 	  dimension u(mx,my+1),v(mx+1,my),p(mx+1,my+1),un(mx,my+1)
      double precision miu
	  miu=1.0/re
	  aw=miu+max(0.5*(u(i-1,j)+u(i,j))*dy,0.0)
      ae=miu+max(0.0,-0.5*(u(i,j)+u(i+1,j))*dy)
      as=miu+max(0.5*(v(i,j-1)+v(i+1,j-1))*dx,0.0)
	  an=miu+max(0.0,-0.5*(v(i,j)+v(i+1,j))*dx)
	  df=0.5*(u(i+1,j)-u(i-1,j))*dy+0.5*(v(i,j)+v(i+1,j)-v(i,j-1)-v(i+1,j-1))*dx
	  ap=aw+ae+as+an+df
	  un(i,j)=u(i,j)+dt/dx/dy*(-ap*u(i,j)+aw*u(i-1,j)+ae*u(i+1,j)  &
+as*u(i,j-1)+an*u(i,j+1))-dt*(p(i+1,j)-p(i,j))/dx
	  Endsubroutine
!
!
      subroutine upbound_v(u,v,p,vn,mx,my,dx,dy,dt,re,i,j)
!��һ��ӭ������ɢ��ʽ�õ��ڲ�߽�vnֵ
	  implicit double precision(a-h,o-z)
 	  dimension u(mx,my+1),v(mx+1,my),p(mx+1,my+1),vn(mx+1,my)
      double precision miu
	  miu=1.0/re
	  aw=miu+max(0.5*(u(i-1,j)+u(i-1,j+1))*dy,0.0)
	  ae=miu+max(0.0,-0.5*(u(i,j)+u(i,j+1))*dy)
	  as=miu+max(0.5*(v(i,j-1)+v(i,j))*dx,0.0)
	  an=miu+max(0.0,-0.5*(v(i,j)+v(i,j+1))*dx)
	  df=0.5*(u(i,j)+u(i,j+1)-u(i-1,j)-u(i-1,j+1))*dy+0.5*(v(i,j+1)-v(i,j-1))*dx
 	  ap=aw+ae+as+an+df
	  vn(i,j)=v(i,j)+dt/dx/dy*(-ap*v(i,j)+aw*v(i-1,j)+ae*v(i+1,j)  &
+as*v(i,j-1)+an*v(i,j+1))-dt*(p(i,j+1)-p(i,j))/dy
   	  Endsubroutine
!
!
      subroutine calp(p,un,vn,pn,mx,my,dx,dy,dt,c2)
!����un��vn���ѹǿpn
      implicit double precision(a-h,o-z)
      dimension p(mx+1,my+1),un(mx,my+1),vn(mx+1,my),pn(mx+1,my+1)
      do i=2,mx
         do j=2,my
	        pn(i,j)=p(i,j)-dt*c2/dx*(un(i,j)-un(i-1,j)+vn(i,j)-vn(i,j-1))
	     enddo
      enddo
      do i=2,mx
         pn(i,1)=pn(i,2)
	     pn(i,my+1)=pn(i,my)
      enddo
      do j=1,my+1
         pn(1,j)=pn(2,j)
	     pn(mx+1,j)=pn(mx,j)
      enddo
      endsubroutine
!
!
      subroutine check(u,v,p,un,vn,pn,mx,my,dx,dy,dt,c2,err)
!У��������Ϣ���õ������ж�׼��err��ͬʱ����u��v��p
	  implicit double precision(a-h,o-z)
	  double precision u(mx,my+1),v(mx+1,my),p(mx+1,my+1)
		  double precision  un(mx,my+1),vn(mx+1,my),pn(mx+1,my+1)
	  do i=1,mx
	     do j=1,my+1
		    abc=abs(un(i,j)-u(i,j))/dt
			if(abc.gt.err) err=abc
		 	u(i,j)=un(i,j)
		 enddo
	  enddo
	  do i=1,mx+1
	     do j=1,my
		    abc=abs(vn(i,j)-v(i,j))/dt
			if(abc.gt.err) err=abc
		    v(i,j)=vn(i,j)
		 enddo
	  enddo
	  do i=1,mx+1
	     do j=1,my+1
		    abc=abs(pn(i,j)-p(i,j))/c2/dt
			if(abc.gt.err) err=abc
		    p(i,j)=pn(i,j)
		 enddo
	  enddo
	  endsubroutine
!
!
   subroutine output(u,v,p,mx,my,dx,dy)
!������
   implicit double precision(a-h,o-z)
  double precision  u(mx,my+1),v(mx+1,my),p(mx+1,my+1),un(mx,my+1)
  double precision  uc(mx,my),vc(mx,my),pc(mx,my),x(mx),y(my)
   do i=1,mx
      x(i)=(i-1)*dx
   enddo
   do j=1,my
      y(j)=(j-1)*dy
   enddo
   do i=1,mx
      do j=1,my
	     uc(i,j)=0.5*(u(i,j)+u(i,j+1))
		 vc(i,j)=0.5*(v(i,j)+v(i+1,j))
	     pc(i,j)=0.25*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))
	  enddo
   enddo
   open(1,file='out.plt')
   write(1,*) 'TITLE = "result"'
   write(1,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
   write(1,*) 'ZONE I=',mx, 'J=',my, 'F=POINT'
   write(1,10) ((x(i),y(j),uc(i,j),vc(i,j),pc(i,j),i=1,mx),j=1,my)
   close(1)
10 format(5(e15.8,5x))
   Endsubroutine
!
!
