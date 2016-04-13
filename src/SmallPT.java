import java.io.*;
import java.util.concurrent.*;
public class SmallPT {
  class Vec {
	final double x, y, z;
	Vec(double ...e) {x = e[0];y = e[1];z = e[2];}
	Vec add(Vec v) { return new Vec(v.x+x,v.y+ y,v.z+z); }
	Vec sub(Vec v) { return new Vec(x-v.x,y-v.y,z-v.z); }
	Vec mul(Vec v) { return new Vec(x*v.x,y*v.y,z*v.z); }
	Vec mul(double d) { return new Vec(d*x,d*y,d*z); }
	Vec n() { return this.mul(1.0/Math.sqrt(x*x+y*y+z*z)); }
	Vec c(Vec v) { return new Vec(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x);}
	double d(Vec v) {return x*v.x+y*v.y+z*v.z;}
  }
  class Ray { 
    Sphere s;
	double dis=Double.MAX_VALUE,t=Double.MAX_VALUE;
	Vec o,d; 
	Ray(Vec o_, Vec d_) {o=o_;d=d_;}
	void intersect() {
	  for (int i=0;i<spheres.length;++i)
		if ((dis=spheres[i].intersect(this))>0&&dis<t) {t=dis;s=spheres[i];}
	}
  }
  enum Refl {DIFF, SPEC, REFR}
  int w=1024,h=768,tot=w*h,p=Runtime.getRuntime().availableProcessors(),samp;
  final Ray cam=new Ray(new Vec(50,52,295.6),new Vec(0,-0.042612,-1));
  final Vec cx=new Vec(w*0.5135/h,0,0),cy =(cx.c(cam.d)).n().mul(0.5135),
		  r=new Vec(.75,.25,.25),bl=new Vec(0,0,0), b=new Vec(.25,.25,.75),
		  gr=new Vec(.75,.75,.75),white=new Vec(0.999,0.999,0.999);
  final Vec[] im = new Vec[w*h];
  public SmallPT(String[] args) throws Exception {
	samp=args.length>0?Math.max(Integer.parseInt(args[0])/4,1):16/4;
	final ExecutorService service = Executors.newFixedThreadPool(p);
	for (int y=0;y<h;++y) service.submit(new Thread(""+y) {
	  public void run() {
		int y=Integer.parseInt(getName());
	    java.util.Random r = new java.util.Random(y*w);
	    for(int x=0;x<w;++x) {
	      Vec t=new Vec(0,0,0);
		  for (int ss=0;ss<4*samp;++ss) {
	        double u=2*r.nextDouble(),dx=u<1?Math.sqrt(u)-1:1-Math.sqrt(2-u);
		    double v=2*r.nextDouble(),dy=v<1?Math.sqrt(v)-1:1-Math.sqrt(2-v);
		    Vec d=cx.mul(((ss%2+0.5+dx)/2.0+x)/w-0.5).add(
			      cy.mul((((ss/2)%2+0.5+dy)/2.0+y)/h-0.5)).add(cam.d);
		    t=t.add(radiance(new Ray(cam.o.add(d.mul(130)),d.n()),0,r));
		  }
		  im[x+w*(h-y-1)]=t.mul(1.0/(4*samp));
	    }
	  }
	});
	service.shutdown();service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
	FileWriter writer = new FileWriter(new File("out.ppm"));
    writer.write("P3\n"+w+" "+h+"\n255"); 
    for(int i=0;i<w*h;i++)
    	writer.write("\n"+g(im[i].x)+" "+g(im[i].y)+" "+g(im[i].z)+"");
  }
  public static void main(String[] a) throws Exception { new SmallPT(a); }
  Sphere[] spheres = {
	new Sphere(1e5,Refl.DIFF,new Vec(1e5+1,40.8,81.6),bl,r),
	new Sphere(1e5,Refl.DIFF,new Vec(-1e5+ 99,40.8,81.6),bl,b),
	new Sphere(1e5,Refl.DIFF,new Vec(50,40.8,1e5),bl,gr),
	new Sphere(1e5,Refl.DIFF,new Vec(50,40.8,-1e5+170),bl,bl),
	new Sphere(1e5,Refl.DIFF,new Vec(50,1e5,81.6),bl,gr),
	new Sphere(1e5,Refl.DIFF,new Vec(50,-1e5+81.6,81.6),bl,gr),
	new Sphere(16.5,Refl.SPEC,new Vec(27,16.5,47),bl,white),
	new Sphere(16.5,Refl.REFR,new Vec(73,16.5,78),bl,white),
	new Sphere(600,Refl.DIFF,new Vec(50,681.33,81.6),new Vec(12,12,12),bl) 
  };
  int g(double x) {return (int)(Math.pow(x<0?0:(x>1?1:x),1.0/2.2)*255+0.5);}
  Vec radiance(Ray r, int depth,java.util.Random ra) {
	r.intersect();
	if (r.s==null) return new Vec(0,0,0);
	Vec x=r.o.add(r.d.mul(r.t)),n=(x.sub(r.s.p)).n(),nl=n.d(r.d)<0?n:n.mul(-1);
	Vec f=r.s.c;
	double p = f.x>f.y&&f.x>f.z?f.x:f.y>f.z?f.z:f.z;
	if (++depth>20) return new Vec(0,0,0); else if (depth>5)
	  if (ra.nextDouble()<p) f=f.mul(1.0/p); else return r.s.e;
	if (r.s.re == Refl.DIFF) {
	  double r1=2*Math.PI*ra.nextDouble(),r2=ra.nextDouble();
	  double r2s=Math.sqrt(r2);
	  Vec w=nl,u=(Math.abs(w.x)>0.1?new Vec(0,1,0):new Vec(1,0,0)).c(w).n();
	  Vec d=((u.mul(Math.cos(r1)*r2s)).add(w.c(u).mul(Math.sin(r1)*r2s)).
			  add(w.mul(Math.sqrt(1-r2)))).n();
	  return r.s.e.add(f.mul(radiance(new Ray(x, d), depth,ra)));
	} else if (r.s.re == Refl.SPEC)
	  return r.s.e.add(f.mul(radiance(
					new Ray(x, r.d.sub(n.mul(2.0 * n.d(r.d)))), depth,ra)));
    Ray reflRay = new Ray(x, r.d.sub(n.mul(2 * n.d(r.d))));
    int in = n.d(nl)>0?1:-1;
    double nc=1,nt=1.5,nnt=in>0?nc/nt:nt/nc,ddn=r.d.d(nl), cos2t;
    if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)
      return r.s.e.add(f.mul(radiance(reflRay, depth,ra)));
    Vec tDir=(r.d.mul(nnt).sub(n.mul(in*(nnt*ddn+Math.sqrt(cos2t))))).n();
    double a=nt-nc,b=nt+nc,R0=a*a/(b*b),c=1-(in>0?-ddn:tDir.d(n));
    double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=0.25+0.5*Re,RP=Re/P,TP=Tr/(1-P);
    if (depth > 2)
      if (ra.nextDouble() < P)
	    return r.s.e.add(f.mul(radiance(reflRay, depth,ra).mul(RP))); else
	    return r.s.e.add(f.mul(radiance(new Ray(x,tDir),depth,ra).mul(TP)));
  	else
	    return r.s.e.add(f.mul(radiance(reflRay, depth,ra)).mul(Re)
				     .add(radiance(new Ray(x, tDir), depth,ra).mul(Tr)));
  }
  class Sphere {
	double rad;
	Vec p,e,c;
	Refl re;
	Sphere(double r, Refl re_,Vec ...v) {rad=r;p=v[0];e=v[1];c=v[2];re=re_;}
	double intersect(Ray r) {
	  Vec op=p.sub(r.o);
	  double t,eps=1e-4,b=op.d(r.d),det=b*b-op.d(op)+rad*rad;
	  return det<0?0:(t=b-(det=Math.sqrt(det)))>eps?t:((t=b+det)>eps?t:0);
	}
  }
}