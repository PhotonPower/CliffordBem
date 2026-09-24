from mfs import *
for (nt,nphi,ns) in [(16,32,50),(20,40,70),(24,48,100),(30,60,140)]:
    s=Sphere(nt,nphi,ns)
    print(nt,nphi,ns, ["%.5f"%min_angle(s,kd,om,e1,1,1.0,1) for kd in ['maxwell','dirac'] for (om,e1) in [(1.0,2.25),(2.0,-11+1.2j)]], flush=True)
