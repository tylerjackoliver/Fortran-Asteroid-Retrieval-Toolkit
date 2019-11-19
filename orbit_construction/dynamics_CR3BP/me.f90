! Set of scripts to single-shoot some orbits in about a lagrangian point 
! L_i 
! 
! TODO: HOW THE FUCK DID I DO THE LINEARISATION?! 
! Variable declarations 
real(MK) :: amplitude1  !<  
real(MK) :: amplitude2  !<  
integer :: initialConds  !<  
real(MK) :: lDist  !<  
real(MK) :: lLoc  !<  
real(MK) :: lPoint  !<  
real(MK) :: massParameter  !<  
integer :: numberOfFamilies  !<  
real(MK) :: OPTIONS  !<  
real(MK) :: T  !<  
real(MK) :: t  !<  
real(MK) :: tb  !<  
real(MK) :: tf  !<  
real(MK) :: x  !<  
real(MK) :: x0  !<  
! 

massParameter = 0.0121 ! System mass parameter
lPoint = 1 ! Desired L-point
numberOfFamilies = 10 ! Number of families to plot
lDist = (massParameter/3)**(1/3) ! Approximate distance to Li from secondary (Ross et. al, 2008)
amplitude1 = 1e-02*lDist
amplitude2 = 2*amplitude1 ! Get the first two sets of amplitudes for continuation procedure
call ( ) [x0, T] = generateInitialGuess(massParameter, lPoint, amplitude1, ... ! Generate initial guess using (Richardson, 1980)
amplitude2, numberOfFamilies)
OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14)
! Set default integrator options
! for k = 1:numberOfFamilies, 
! initialConds = x0(k,:);						 ! Extract initial conditions from the IC array 
! tb = 0;								 ! Set the integration time to begin at x0 
! tf = T(k)/2;								 ! Integrate half of the orbit 
! [x,t] = integrator(initialConds, tb, tf, massParameter, OPTIONS); 
! plot(x(:,1),x(:,2),'k',x(:,1),-x(:,2),'b');				 ! Plot result 
! hold on 
! axis equal 
! [lLoc] = lPointFinder(massParameter, lPoint);				 ! Append lagrangian point to plot -- find location 
! plot(lLoc(1), lLoc(2), '-ro'); 
! end 
initialConds = x0(1,:)
! figure(); 
tf = T(1)
tb = 0
call integrator( initialConds, tb, tf, massParameter, OPTIONS , x,t) 
plot(x(:,1),x(:,2),'k')
axis equal
call lPointFinder( massParameter, lPoint , lLoc) 
plot(lLoc(1), lLoc(2), 'MarkerSize', 10)
!> 
subroutine generateInitialGuess(x0, T)
implicit none 
! Arguments declarations 
real(MK), dimension(2), intent(out) :: T  !< ! m2f:check dim(numb 
real(MK), dimension(:,:), intent(out) :: x0  !< ! m2f:check dim(numb 
! Variable declarations 
real(MK) :: dimensions  !<  
real(MK) :: dt  !<  
real(MK) :: dx  !<  
real(MK) :: dydot_product  !<  
real(MK) :: familyNumber  !<  
real(MK) :: t1  !<  
real(MK) :: t2  !<  
real(MK) :: tempState  !<  
real(MK) :: tempTime  !<  
real(MK) :: x01  !<  
real(MK) :: x02  !<  
real(MK) :: x0Guess  !<  
real(MK) :: x0Guess1  !<  
real(MK) :: x0Guess2  !<  
! 
amplitude2, numberOfFamilies)
dt = -1e-02 ! Guess at change in period between orbits
dimensions = 4 ! Number of dimensions we're dealing with (x(dot_product), y(dot_product))
!m2f: x0 = zeros(numberOfFamilies, dimensions) ! Initialise matrices for states and periods
if (allocated(x0)) deallocate(x0)
allocate(x0(numberOfFamilies, dimensions))
x0 = 0.0_MK 

!m2f: T = zeros(numberOfFamilies, 1)
if (allocated(T)) deallocate(T)
allocate(T(numberOfFamilies, 1))
T = 0.0_MK 


x0Guess1 = initialGuessGenerator(massParameter, lPoint, amplitude1) ! Send these off for an initial guess(timate)
x0Guess2 = initialGuessGenerator(massParameter, lPoint, amplitude2)

! The following methodology is shamelessly copied from Ross, 2008: use 
! dx/dy from two different orbits in the family to seed the next set 

familyNumber = 1
!m2f: disp(sprintf('Getting orbit family:  write(s,*) ' '  
call differentialCorrection( x0Guess1, massParameter , x01,t1) ! Send off the initial conditions to get refined guess

familyNumber = 2
!m2f: disp(sprintf('Getting orbit family:  write(s,*) ' '  
call differentialCorrection( x0Guess2, massParameter , x02,t2) 

x0(1:2, 1:dimensions) = x01(:)'; x02(:)';! Append the 'good' initial conditions to the matrix of init. conds.
T(1:2) = [ 2*t1,2*t2 ] ! Transpose to get 4x1 -> 1x4

do familyNumber=3,numberOfFamilies ! Gotta catch 'em all! 
!m2f: disp(sprintf('Getting orbit family:  write(s,*) ' '  
dx = x0(familyNumber-1, 1) - x0(familyNumber-2, 1) ! Compute the deltas between the last two calc'd orbits
dydot_product = x0(familyNumber-1, 4) - x0(familyNumber-2, 4)
dt = T(familyNumber-1) - T(familyNumber-2)
x0Guess = [(x0(familyNumber-1,1)+dx) 0 0 (x0(familyNumber-1,&
4)+dydot_product)]
call ( ) [tempState, tempTime] = differentialCorrection(x0Guess, ...
massParameter)
x0(familyNumber, 1:dimensions) = tempState
T(familyNumber) = 2*tempTime ! we only integrate 1/2, so x2!
end do 
end subroutine generateInitialGuess 


!> 
subroutine pointDistance(lPoint, massParameter, gamma)
implicit none 
! Arguments declarations 
real(MK), intent(out) :: gamma  !<  
real(MK), intent(in) :: lPoint  !<  
real(MK), intent(in) :: massParameter  !<  
! Variable declarations 
real(MK), dimension(6) :: l1Poly  !<  
real(MK) :: l1Roots  !<  
real(MK), dimension(6) :: l2Poly  !<  
real(MK) :: l2Roots  !<  
real(MK), dimension(6) :: l3Poly  !<  
real(MK) :: l3Roots  !<  
real(MK) :: m1  !<  
real(MK) :: m2  !<  
! 
m1 = massParameter ! Dependecy function for locating position of lPoint
m2 = 1-massParameter

! Solve the quintic polynomial found in (Szebhely, 1967) 
l1Poly = [ 1,(m1-3),(3-2*m1),-m1,2*m1,-m1 ] 
l2Poly = [ 1,(3-m1),(3-2*m1),-m1,-2*m1,-m1 ] 
l3Poly = [ 1,(2+m1),(1+2*m1),-m2,-2*m2,-m2 ] 

l1Roots = roots(l1Poly)
l2Roots = roots(l2Poly)
l3Roots = roots(l3Poly)

! Get rid of the complex pairs 
do k=1,5 
if (isreal(l1Roots(k)) gammas(1)=l1Roots(k)) then 
end if 
if (isreal(l2Roots(k)) gammas(2)=l2Roots(k)) then 
end if 
if (isreal(l3Roots(k)) gammas(3)=l3Roots(k)) then 
end if 
end do 
gamma = gammas(lPoint)
end subroutine pointDistance 

!> 
subroutine lPointFinder(massParameter, lPoint, lPointPos)
implicit none 
! Arguments declarations 
real(MK), intent(in) :: lPoint  !<  
real(MK), intent(out) :: lPointPos  !<  
real(MK), intent(in) :: massParameter  !<  
! Variable declarations 
real(MK), dimension(:,:) :: lPoints  !< ! m2f:check dim(5, 2 
real(MK) :: m1  !<  
real(MK) :: m2  !<  
real(MK) :: tempLPoints  !<  
! 
m1 = 1-massParameter
m2 = massParameter
!m2f: lPoints = zeros(5, 2) ! Set up position matrix in planar case
if (allocated(lPoints)) deallocate(lPoints)
allocate(lPoints(5, 2))
lPoints = 0.0_MK 

tempLPoints = [(m1 - pointDistance(1, massParameter)) & ! Get the position of the three collienar points
(m1 + pointDistance(2, massParameter)) (-m2&
-pointDistance(3, massParameter))]
lPoints(1, 1) = tempLPoints(1) ! Three collinear points are at y=0, so no need to assign from zeros matrix
lPoints(2, 1) = tempLPoints(2)
lPoints(3, 1) = tempLPoints(3)
lPoints(4, 1) = 0.5-massParameter ! X-coord of equilateral solutions
lPoints(5, 1) = 0.5-massParameter
lPoints(4, 2) = 0.5*sqrt(3) ! Y-coord of equilateral solutions
lPoints(5, 2) = -0.5*sqrt(3)

lPointPos = lPoints(lPoint, :) ! Return desired position
end subroutine lPointFinder 

!> 
subroutine jacobianMatrix(lPos, massParameter, jacobian)
implicit none 
! Arguments declarations 
real(MK), dimension(4,4), intent(out) :: jacobian  !<  
real(MK), intent(in) :: lPos  !<  
real(MK), intent(in) :: massParameter  !<  
! Variable declarations 
real(MK) :: m1  !<  
real(MK) :: m2  !<  
real(MK) :: r2  !<  
real(MK) :: R2  !<  
real(MK) :: r3  !<  
real(MK) :: R3  !<  
real(MK) :: r5  !<  
real(MK) :: R5  !<  
real(MK) :: uxx  !<  
real(MK) :: uxy  !<  
real(MK) :: uyy  !<  
real(MK) :: x  !<  
! 
m1 = 1-massParameter
x = lPos
m2 = massParameter
r2= (x(1)+m2)**2+x(2)**2 ! r: distance to m1, LARGER MASS
R2= (x(1)-m1)**2+x(2)**2 ! R: distance to m2, smaller mass
! Calculate powers to ease up our life a little 
r3= r2**1.5
r5= r2**2.5
R3= R2**1.5
R5= R2**2.5
uxx = -1+(m1/r3)*(1-(3*(x(1)+m2)**2/r2))+(m2/R3)*(1-(3*(x(1)-m1)**2 & ! Compute the derivatives of the potential U
/R2))
uyy = -1+(m1/r3)*(1-(3*x(2)**2/r2))+(m2/R3)*(1-(3* x(2)**2/R2))
uxy = -(m1/r5)*3*x(2)*(x(1)+m2)-(m2/R5)*3*x(2)*(x(1)-m1)
jacobian = reshape([ 0,0,1,0,0,0,0,1,-uxx,-uyy,0,2,-uxy,-uyy,-2,0 ] , [ 4 , 4] )! Arrange into restricted jacobian for planar case
end subroutine jacobianMatrix 

!> 
subroutine getEigenvalues(matrix, eValStable, eValUnstable, eValCenter)
implicit none 
! Arguments declarations 
real(MK), dimension(1), intent(out) :: eValCenter  !<  
real(MK), dimension(1), intent(out) :: eValStable  !<  
real(MK), dimension(1), intent(out) :: eValUnstable  !<  
real(MK), intent(in) :: matrix  !<  
! Variable declarations 
real(MK) :: centerCounter  !<  
real(MK) :: directions  !<  
real(MK) :: height  !<  
character(len=*) :: stableCounter  !<  
real(MK) :: unstableCounter  !<  
real(MK), dimension(:) :: values  !<  
real(MK) :: width  !<  
! 
eValStable = [ ] 
eValUnstable = [ ] 
eValCenter = [ ] 
call shape( matrix , height,width) ! We need to analyse each bit of the returned matrix, so grab the dimensions
call eig( matrix , values,directions) ! Compute matrix eigenvals/vecs (not bothered about vec atm)
stableCounter = 0 ! We will categorise the nature of the eigenvalue, so initialise some counters
unstableCounter = 0
centerCounter = 0
directions = cleanUpMatrix(directions)
do i=1,height 
if (real(directions(i, i)) < 0) then ! Measure stability using sign
stableCounter = stableCounter + 1 ! Increment counter and add to array
eValStable = directions(i, i)
else if (real(directions(i, i)) > 0) then 
unstableCounter = unstableCounter + 1
eValUnstable = directions(i, i)
else
centerCounter = centerCounter + 1
eValCenter = directions(i, i)
end if 
end do 
end subroutine getEigenvalues 

!> 
subroutine lPointEigenvalues(eValStable, eValUnstable, eValCenter)
implicit none 
! Arguments declarations 
real(MK), intent(out) :: eValCenter  !<  
real(MK), intent(out) :: eValStable  !<  
real(MK), intent(out) :: eValUnstable  !<  
! Variable declarations 
real(MK) :: jacobian  !<  
! 
massParameter) ! Catch-all function to generate the eigvenVals about a point
jacobian = jacobianMatrix(lPos, massParameter)
call getEigenvalues( jacobian , eValStable,eValUnstable,eValCenter) 
end subroutine lPointEigenvalues 

!> 
subroutine initialGuessGenerator(x0Guess, tGuess)
implicit none 
! Arguments declarations 
real(MK), intent(out) :: tGuess  !<  
real(MK), dimension(:), intent(out) :: x0Guess  !< ! m2f:check dim(4, 1 
! Variable declarations 
real(MK) :: c2  !<  
real(MK) :: centerEigenvalueFrequency  !<  
real(MK) :: eValCenter  !<  
real(MK) :: eValStable  !<  
real(MK) :: eValUnstable  !<  
real(MK) :: gamma  !<  
real(MK) :: k  !<  
real(MK) :: l  !<  
real(MK), dimension(4) :: lPos  !<  
! 
lPoint, amplitude) ! Generate an initial guess using (Richardson, 1980
Halo Orbit Formulation for the ISEE-3 mission)
lPos = [ lPointFinder(massParameter,lPoint),0,0 ] ! Set up position of collinear lagrange point
call ( ) [eValStable, eValUnstable, eValCenter] = lPointEigenvalues(lPos,... ! Get the eigenvalues for seeding initial guess from center eval
massParameter)
centerEigenvalueFrequency = abs(imag(eValCenter(1))) ! Imaginary part of the center eigenvalue is the frequency
gamma = pointDistance(lPoint, massParameter) ! Get distance to L_i (Szebhely, 1967)
l = centerEigenvalueFrequency(1) ! We want the frequency of the center eigenvalue, to avoid perturbing direction of propagation
c2 = gamma**(-3)*(massParameter + (-1)**2*(1-massParameter)*&
(gamma/(1-gamma))**(2+1))
k = (l**2+1+2*c2)/(2*l)
!m2f: x0Guess = zeros(4, 1) ! Initialise state matrix
if (allocated(x0Guess)) deallocate(x0Guess)
allocate(x0Guess(4, 1))
x0Guess = 0.0_MK 

x0Guess(1) = lPos(1) - amplitude ! Set up x guess
x0Guess(4) = lPos(4) + amplitude*k*l ! Set up y-dot_product guess
tGuess = 2*pi/l ! Guess the time
end subroutine initialGuessGenerator 

!> 
subroutine differentialCorrection(x0, massParameter, time1)
implicit none 
! Arguments declarations 
real(MK), intent(in) :: massParameter  !<  
real(MK), intent(out) :: time1  !<  
real(MK), dimension(:), intent(inout) :: x0  !< !m2f: check dim(:) 
! Variable declarations 
real(MK) :: AbsTol  !<  
real(MK) :: attempt  !<  
real(MK) :: dampFactor  !<  
real(MK) :: m  !<  
real(MK) :: m1  !<  
real(MK) :: m2  !<  
real(MK) :: MAXdxdot_product1  !<  
real(MK) :: maximumAttempts  !<  
real(MK) :: model  !<  
real(MK) :: options  !<  
real(MK) :: phi_t1  !<  
real(MK) :: r3  !<  
real(MK) :: R3  !<  
real(MK) :: RelTol  !<  
real(MK) :: temp1  !<  
real(MK) :: temp2  !<  
real(MK) :: temp3  !<  
real(MK), dimension(2) :: time  !<  
real(MK) :: u_x  !<  
real(MK) :: x  !<  
real(MK) :: x1  !<  
real(MK) :: xAccel  !<  
real(MK) :: xdot_product1  !<  
real(MK) :: xx  !<  
real(MK) :: y1  !<  
real(MK), dimension(:) :: yAccelDelta  !< !m2f: check dim(:) 
real(MK) :: ydot_product1  !<  
! 
MAXdxdot_product1 = 1e-010
RelTol = 3.e-14
AbsTol = 1.e-16 ! Integration options
maximumAttempts = 100
xdot_product1 = 1
attempt = 0

do while (abs(xdot_product1) > MAXdxdot_product1) 
if (attempt > maximumAttempts) then ! Computation handling
disp('Oh shit, too many iterations!')
exit
end if 
model = 'cr3bp' ! Point MATLAB to the governing equations
time = [ 0,10 ] ! Arbritrary time span for the integration
options = odeset('RelTol', RelTol, 'AbsTol', AbsTol, 'Events', 'on')
call ode113( model, time, x0, options , temp1,temp2,time1,xx,temp3) ! Integrate until x-axis crossing
x1 = xx(1) ! Extract state at x-axis crossing
y1 = xx(2)
xdot_product1 = xx(3)
ydot_product1 = xx(4)
if (time1(1) < 1e-01) then 
time1 = time1(2)
end if 
options = odeset('RelTol', RelTol, 'AbsTol', AbsTol)
call ( ) [x, temp1, phi_t1, temp2] = stateTransitionMatrix(x0, time1, ... ! Propagate the STM over the orbit
massParameter, options)
attempt = attempt + 1
!m2f: disp(fprintf('Differential correction: hoping this works number  write(unitnumber,*) ' ' 
attempt))
plot(x(:,1),x(:,2),'r-')
hold on
m = size(x)
plot(x(1,1),x(1,2))
plot(x(m,1), x(1,2))
pause(0.01)

! Compute change in state for next attempt at correction 
m1 = 1-massParameter ! normalised masses of primaries in cr3bp
m2 = massParameter
r3 = ((x1+m2)**2+y1**2)**1.5 ! r: distance to m1
R3 = ((x1-m1)**2+y1**2)**1.5 ! R: distance to m2
u_x = -x1+m1*(x1+m2)/r3 + m2*(x1-m1)/R3 ! x-deriv. of the potential gunction
xAccel = 2*ydot_product1-u_x ! Compute the change in xAccel
yAccelDelta(attempt) = &
(1/(phi_t1(3,4)-phi_t1(2,4)*(1/ydot_product1)*xAccel))*xdot_product1 ! Compute the change in yAccel ! Apply the changes
if (massParameter > 1e-03) then 
!dampFactor = 1-0.5^attempt; 
dampFactor = 1
else
dampFactor = 1
end if 
x0(4) = x0(4) - dampFactor*yAccelDelta(attempt)
end do 
disp('FUCK YEA, THAT WORKED')
end subroutine differentialCorrection 

!> 
subroutine stateTransitionMatrix(x, time, phi, PHI)
implicit none 
! Arguments declarations 
real(MK), intent(out) :: phi  !<  
real(MK), intent(out) :: PHI  !<  
real(MK), intent(out) :: time  !<  
real(MK), intent(out) :: x  !<  
! Variable declarations 
real(MK) :: dimension  !<  
real(MK) :: fixed_step  !<  
real(MK) :: model  !<  
real(MK) :: OPTIONS  !<  
real(MK) :: phi  !<  
real(MK), dimension(:) :: PHI_0  !< !m2f: check dim(:)!m 
real(MK), dimension(2) :: time1  !<  
real(MK) :: time2  !<  
! 
massParameter, options, fixedStepFlag)
dimension = 4 ! Dimension of the problem; see above
model = 'varEqs3BP' ! Point MATLAB to the variational equations model
if (nargin < 5) then ! Argument control for fixedStepFlag
fixed_step = 0
if (nargin < 4) then 
OPTIONS = odeset('RelTol',3e-14,'AbsTol',1e-14) ! integration options, if not specified
end if 
end if 

if (fixed_step == 0) then 
time1 = [ 0,time ] 
else
time1 = 0:time/(fixed_step-1):time
end if 
PHI_0(1:dimension**2) = reshape(eye(dimension), dimension**2, 1) ! Initialise stm to (16x1) identity matrix
PHI_0(1+dimension**2:dimension+dimension**2) = x0 ! Initialise the size(XXX,XXX) of PHI to be the state
call ode113( model, time1, PHI_0, options , time2,PHI) ! Integrate the STM over the full orbit
x = PHI(:,1+dimension**2:dimension+dimension**2)
! trajectory from time 0 to tf ! Re-grab trajectory
phi = reshape(PHI(size(time2), 1:dimension**2), dimension, dimension) ! Grab the actual STM
end subroutine stateTransitionMatrix 


!> 
subroutine integrator(x0, tb, tf, massParameter, options, x, t)
implicit none 
! Arguments declarations 
real(MK), intent(in) :: massParameter  !<  
real(MK), intent(in) :: options  !<  
real(MK), dimension(2), intent(out) :: t  !<  
real(MK), intent(in) :: tb  !<  
real(MK), intent(in) :: tf  !<  
real(MK), dimension(:), intent(out) :: x  !<  
real(MK), intent(in) :: x0  !<  
! Variable declarations 
real(MK) :: model  !<  
real(MK), dimension(2) :: timetb  !<  
real(MK), dimension(2) :: timetf  !<  
real(MK) :: tt1  !<  
real(MK) :: tt2  !<  
real(MK) :: xx1  !<  
real(MK) :: xx2  !<  
! 
model = 'cr3bp'
options = odeset('RelTol',3e-10,'AbsTol',1e-12, 'Events', 'off')
if (nargin < 5) then ! Option handling
options = odeset('RelTol',3e-10,'AbsTol',1e-12, 'Events', 'off')
end if 
timetb = [ -1e-012,tb ] 
timetf = [ 0,tf ] 
x = [ ] ! Initialise output array
if (tb == 0) then 
call ode113( model, timetf, x0, options , t,x) 
else if (tf == 0) then 
call ode113( model, timetb, x0, options , t,x) 
else
call ode113( model, timetb, x0, options , tt1,xx1) 
call ode113( model, timetb, x0, options , tt2,xx2) 
x = [ flipud(xx1),xx2 ] 
t = [ flipud(tt1),tt2 ] 
end if 
end subroutine integrator 


!> 
subroutine cleanUpMatrix(A)
implicit none 
! Arguments declarations 
real(MK), dimension(:,:), intent(inout) :: A  !< !m2f: check dim(:,:) 
! Variable declarations 
real(MK) :: tolerance  !<  
! 
tolerance = 1e-014
do k=1,size(A) 
do l=1,size(A) 
if (abs(real(A(k,1))) < tolerance) then 
A(k,l) = i*imag(A(k,l))
end if 
if (abs(imag(A(k,l))) < tolerance) then 
A(k,l) = real(A(k,l))
end if 
end do 
end do 
end subroutine cleanUpMatrix 

!> 
subroutine manifoldFinder(massParameter, x0, t0, dimensions, N, manifoldIC)
implicit none 
! Arguments declarations 
real(MK), intent(in) :: dimensions  !<  
real(MK), dimension(4), intent(out) :: manifoldIC  !<  
real(MK), intent(in) :: massParameter  !<  
real(MK), intent(in) :: N  !<  
real(MK), intent(in) :: t0  !<  
real(MK), intent(in) :: x0  !<  
! Variable declarations 
real(MK) :: .not.(  !<  
real(MK) :: a  !<  
integer :: i1  !<  
integer :: i2  !<  
integer, dimension(2) :: initialConds  !<  
real(MK), dimension(1) :: M  !<  
real(MK) :: monodromyEigVec  !<  
real(MK) :: OPTIONS  !<  
real(MK), dimension(1) :: Phi  !<  
real(MK) :: real  !<  
real(MK) :: S  !<  
character(len=*) :: stable  !<  
character(len=*) :: stm  !<  
real(MK) :: U  !<  
real(MK) :: unstable  !<  
real(MK), dimension(:) :: val  !<  
real(MK) :: X  !<  
real(MK) :: x  !<  
real(MK), dimension(:,:) :: xs_in  !< ! m2f:check dim(dime 
real(MK), dimension(:,:) :: xs_out  !< !m2f: check dim(:,:) 
real(MK), dimension(:,:) :: xu_in  !< !m2f: check dim(:,:) 
real(MK), dimension(:,:) :: xu_out  !< !m2f: check dim(:,:) 
! 
stm = eye(dimensions) ! Initialise STM
initialConds = [ x0,stm ] ! Initialise...initial conditions
OPTIONS = odeset('RelTol', 1e-013, 'AbsTol', 1e-013)
t0 = linspace(0, t0, N)
call ode45( @crtbp, t0, initialConds, OPTIONS, massParameter , .not.(,X) )
M = reshape(X(size(XXX,XXX), 5:size(XXX,XXX)),6,[ ] )! Extract monodromy matrix
call eig( M , monodromyEigVec,val) 
real = (imag(diag(val)) == 0)
call max( abs(val) , .not.(,i1) )
call min( abs(val) , .not.(,i2) )
a = 1:6
stable = monodromyEigVec(:, a(i1))
unstable = monodromyEigVec(:, a(i2))
!m2f: xs_in = zeros(dimensions, N)
if (allocated(xs_in)) deallocate(xs_in)
allocate(xs_in(dimensions, N))
xs_in = 0.0_MK 

xs_out = xs_in
xu_in = xs_in
xu_out = xs_in
do i=1,N 
Phi = reshape(X(i,5:size(XXX,XXX)), [ ] )
x = X(i, 1:4)
S = Phi*stable
U = Phi*unstable
xs_out(:, i) = x(:, i) + S.*10**(-4)
xs_in(:,i) = x(:, i) - S.*10**(-4)
xu_out(:,i) = x(:,i) + U.*10**(-4)
xu_in(:,i) = x(:,i) - U.*10**(-4)
end do 
manifoldIC = [ xs_out,xs_in,xu_out,xu_in ] 
end subroutine manifoldFinder 

end 
