#!/usr/bin/python

'''
Based on The QCP Superposition Method discussed in detail at: http://theobald.brandeis.edu/qcp/

A pure python implementation of the fast rmsd / rotation algorithm using special properties of 
the quaternion representation of rotation. 

The QCP method is the fastest method known for determining the minimum RMSD between two structures 
and for determining the optimal least-squares rotation matrix. Least-squares superposition methods 
find the rotation matrix that minimizes the RMSD. The most common algorithms use an eigendecomposition, 
a singular value decomposition, or an inversion of a small "key" matrix (of rank 3 or 4). The QCP 
method avoids the costly matrix decomposition by taking advantage of special properties of the 
quaternion representation of rotations and of the characteristic polynomial of the 4×4 "key" matrix. 
'''

def fast_calc_rmsd(A, E0, minScore = 0, evalprec = 1e-11):
  C = [0] * 3; 
  Sxx = A[0]; Sxy = A[1]; Sxz = A[2];
  Syx = A[3]; Syy = A[4]; Syz = A[5];
  Szx = A[6]; Szy = A[7]; Szz = A[8];

  Sxx2 = Sxx * Sxx;
  Syy2 = Syy * Syy;
  Szz2 = Szz * Szz;

  Sxy2 = Sxy * Sxy;
  Syz2 = Syz * Syz;
  Sxz2 = Sxz * Sxz;

  Syx2 = Syx * Syx;
  Szy2 = Szy * Szy;
  Szx2 = Szx * Szx;

  SyzSzymSyySzz2 = 2*(Syz*Szy - Syy*Szz);
  Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

  C[2] = -2 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
  C[1] = 8 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

  SxzpSzx = Sxz + Szx;
  SyzpSzy = Syz + Szy;
  SxypSyx = Sxy + Syx;
  SyzmSzy = Syz - Szy;
  SxzmSzx = Sxz - Szx;
  SxymSyx = Sxy - Syx;
  SxxpSyy = Sxx + Syy;
  SxxmSyy = Sxx - Syy;
  Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

  C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2 \
  	+ (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2) \
	+ (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz)) \
	+ (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz)) \
	+ (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz)) \
	+ (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));  
  # Newton-Raphson 
  mxEigenV = E0;
  for i in range(0, 50):
    oldg = mxEigenV;
    x2 = mxEigenV*mxEigenV;
    b = (x2 + C[2])*mxEigenV;
    a = b + C[1];
    delta = ((a*mxEigenV + C[0])/(2*x2*mxEigenV + b + a));
    mxEigenV -= delta;
    if (abs(mxEigenV - oldg) < abs(evalprec*mxEigenV)):
      break;

  if (i == 50):
    print("\nMore than %d iterations needed!\n", i);

  #
  rms = abs(2 * (E0 - mxEigenV));
  rmsd = rms;

  if minScore > 0:
    if rms < minScore:
      return -1; # Don't bother with rotation. 
  
  a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx;
  a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx;
  a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy;
  a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
  a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34;
  a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33;
  a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32;
  q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
  q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
  q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
  q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

  qsqr = q1**2 + q2**2 + q3**2 + q4**2
  
  return (rmsd, (q1, q2, q3, q4))

def inner_product(coords1, coords2, weight = None):
  A = [0] * 9; G1 = 0; G2 = 0
  if weight == None:
    for (x1, y1, z1), (x2, y2, z2) in zip(coords1, coords2):
      G1 += x1**2 + y1**2 + z1**2
      G2 += x2**2 + y2**2 + z2**2
      A[0] +=  (x1 * x2);
      A[1] +=  (x1 * y2);
      A[2] +=  (x1 * z2);
                                    
      A[3] +=  (y1 * x2);
      A[4] +=  (y1 * y2);
      A[5] +=  (y1 * z2);
                                                                        
      A[6] +=  (z1 * x2);
      A[7] +=  (z1 * y2);
      A[8] +=  (z1 * z2);       
  
  return ((G1+G2)/2, A)

def center_coords(coords, weight = None):
  xsum = 0; ysum = 0; zsum = 0; l = len(coords)
  if weight == None:
    for (x, y, z) in coords:
      xsum += x; ysum += y; zsum += z;
    xsum /= l
    ysum /= l
    zsum /= l
  return [(x-xsum,y-ysum,z-zsum) for (x, y, z) in coords]

def align(coords1, coords2, verbose = False):
  centered_coords1 = center_coords(coords1)
  centered_coords2 = center_coords(coords2)
  (E0, A) = inner_product(centered_coords1, centered_coords2)
  (rmsd, quarts) = fast_calc_rmsd(A, E0)
  if verbose:
    print E0
    print A
    print quarts
  ##
  rot = [0] * 9
  qsqr = 0
  for q in quarts:
    qsqr += q**2
  (q1, q2, q3, q4) = quarts
  a2 = q1 * q1 / qsqr;
  x2 = q2 * q2 / qsqr;
  y2 = q3 * q3 / qsqr;
  z2 = q4 * q4 / qsqr;

  xy = q2 * q3 / qsqr;
  az = q1 * q4 / qsqr;
  zx = q4 * q2 / qsqr;
  ay = q1 * q3 / qsqr;
  yz = q3 * q4 / qsqr;
  ax = q1 * q2 / qsqr;

  rot[0] = a2 + x2 - y2 - z2;
  rot[1] = 2 * (xy + az);
  rot[2] = 2 * (zx - ay);
  rot[3] = 2 * (xy - az);
  rot[4] = a2 - x2 + y2 - z2;
  rot[5] = 2 * (yz + ax);
  rot[6] = 2 * (zx + ay);
  rot[7] = 2 * (yz - ax);
  rot[8] = a2 - x2 - y2 + z2;
  return (rmsd, rot, quarts)  
  
def do_test(coords1, coords2, test_no):
  (rmsd, rot, quarts) = align(coords1, coords2, True)
  qsqr = 0
  for q in quarts:
    qsqr += q**2
  norm = qsqr**0.5
  quarts = [q/norm for q in quarts]
  ## Result
  print '-----Test %d-----' %test_no
  print (rmsd/len(coords1))**0.5
  for i in range(0,3):
    print (rot[3*i], rot[3*i+1], rot[3*i+2])
  print quarts
  
if __name__ == '__main__':
  ## Test data
  coords1 = [(-2.803,-15.373,24.556),(0.893,-16.062,25.147),(1.368,-12.371,25.885), \
    (-1.651,-12.153,28.177),(-0.440,-15.218,30.068),(2.551,-13.273,31.372), \
    (0.105,-11.330,33.567)]
  coords2 = [(-14.739,-18.673,15.040),(-12.473,-15.810,16.074),(-14.802,-13.307,14.408), \
    (-17.782,-14.852,16.171),(-16.124,-14.617,19.584),(-15.029,-11.037,18.902), \
    (-18.577,-10.001,17.996)]
  ## Test 1
  do_test(coords1, coords2, 1)
  ## Test 2
  k=10
  coords1 = [(int(c[0]*k),int(c[1]*k),int(c[2]*k)) for c in coords1]
  coords2 = [(int(c[0]*k),int(c[1]*k),int(c[2]*k)) for c in coords2]  
  do_test(coords1, coords2, 2)

      
  