//principle for not reusing momenta on an L^4 lattice

	for(l0=0; l0<L; ++l0)
	{
		for(l1=l0+1; l1<L; ++l1)
		{
			for(l2=l1+1; l2<L; ++l2)
			{
				for(l3=l2+1; l3<L; ++l3)
				{
					All are different - comes 24 times
				}
				{
					l3=l2
					largest 2 are equal - comes 12 times
				}
			}
			{
				l2=l1
				for(l3=l2+1; l3<L; ++l3)
				{
					momenta in the middle are equal - comes 12 times
				}
				{
					l3=l2
					largest 3 momenta are equal - comes 4 times
				}
			}
		}
		{
			l1=l0
			for(l2=l1+1; l2<L; ++l2)
			{
				for(l3=l2+1; l3<L; ++l3)
				{
					smallest 2 are equal - comes 12 times
				}
				{
					l3=l2
					2 x 2 momenta are equal - comes 6 times
				}
			}
			{
				l2=l1
				for(l3=l2+1; l3<L; ++l3)
				{
					smallest 3 momenta are equal - comes 4 times
				}
				{
					//test here for zero mode if neccessary
					l3=l2
					all are equal - comes 1 time
				}
			}
		}
	}

============================================================
============================================================
============================================================

//principle for not reusing momenta on an L^3x2L lattice
//here l0 is used for the larger dimension (the momenta have to be computed correctly of course, I skip that here)

	for(l0=0; l0<L; ++l0)
	{
		//NOTE take the momenta 2pi/2L*(2*l0), 
		for(l1=l0+1; l1<L; ++l1)
		{
			for(l2=l1+1; l2<L; ++l2)
			{
				for(l3=l2+1; l3<L; ++l3)
				{
					All are different - comes 24 times
				}
				{
					l3=l2
					largest 2 are equal - comes 12 times
				}
			}
			{
				l2=l1
				for(l3=l2+1; l3<L; ++l3)
				{
					momenta in the middle are equal - comes 12 times
				}
				{
					l3=l2
					largest 3 momenta are equal - comes 4 times
				}
			}
		}
		{
			l1=l0
			for(l2=l1+1; l2<L; ++l2)
			{
				for(l3=l2+1; l3<L; ++l3)
				{
					smallest 2 are equal - comes 12 times
				}
				{
					l3=l2
					2 x 2 momenta are equal - comes 6 times
				}
			}
			{
				l2=l1
				for(l3=l2+1; l3<L; ++l3)
				{
					smallest 3 momenta are equal - comes 4 times
				}
				{
					//test here for zero mode if neccessary
					l3=l2
					all are equal - comes 1 time
				}
			}
		}
	}
	//now the part where the odd momenta for l0 come (only l1,l2 and l3 cann be commuted)
	for(l0=0; l0<L; ++l0)
	{
		//NOTE take the momenta 2pi/2L*(2*(l0+1)), 
		for(l1=0; l1<L; ++l1)
		{
			for(l2=l1+1; l2<L; ++l2)
			{
				for(l3=l2+1; l3<L; ++l3)
				{
					All three are different - comes 6 times
				}
				{
					l3=l2
					larger 2 are equal - comes 3 times
				}
			}
			{
				l2=l1
				for(l3=l2+1; l3<L; ++l3)
				{
					smaller two are equal - comes 3 times
				}
				{
					l3=l2
					all are equal - comes 1 times
				}
			}
		}
	}

============================================================
============================================================
============================================================

//now the principle for the case, that L/2-n and L/2+n are equivalent
// only consider the L%2=0 case
//the tricky part is, that 0 amd L/2 occur only once
//L^4 lattice

possibilities: p,q,r,s are 1...L/2-1
p<q<r<s

{p,q,r,s}   = 24*16=384             a1
{p,q,r,L/2} = 24*8=192              a2
{0,p,q,r}   = 24*8=192              a3
{0,q,r,L/2} = 24*4=96               a4

{p,p,r,s} = {p,q,q,r} = {p,q,r,r}    = 12*16 = 192     b1(A,B,C)
{p,p,r,L/2} = { p,q,q, L/2}          = 12*8 = 96       b2(A,B)
{0,q,q,r} = {0,q,r,r}                = 12*8 = 96       b3 (A;B)
{0,q,q,L/2}                          = 12*4 =48        b4
{p,q,L/2,L/2}                        = 12*4 = 48       b5
{0,0,r,s}                            = 12* 4 = 48      b6
{0,q,L/2,L/2}                        = 12*2 = 24       b7
{0,0,r,L/2}                          = 12*2 = 24       b8

{p,p,p,s} = {p,q,q,q}             =4*16 =64            c1(A;B)
{p,p,p,L/2}                       = 4*8 = 32           c2
{0,q,q,q}                         = 4*8 = 32           c3
{p,L/2,L/2,L/2}                   = 4*2 = 8            c4
{0,0,0,q}                         = 4*2 = 8            c5
{0,0,0,L/2}                       = 4 *1 =4            c6
{0,L/2,L/2,L/2}                   =  4*1 =4            c7

{p,p,r,r}         = 6 * 16                             d1
{p,p,L/2,L/2)     = 6 * 4                              d2
{0,0,r,r}         = 6 * 4                              d3
{0,0,L/2,L/2}     = 6*1                                d4

{p,p,p,p}        =1*16                                 e1 
{L/2,L/2,L/2,L/2}=1*1                                  e2
{0,0,0,0} = 1*1                                        e3


	for(l0=1; l0<L/2; ++l0)
	{
		for(l1=l0+1; l1<L/2; ++l1)
		{
			for(l2=l1+1; l2<L/2; ++l2)
			{
				for(l3=l2+1; l3<L/2; ++l3)
				{
					All are different - comes 24*16=384 times a1
				}
				{	
					l3=L/2
					all different, one L/2 - comes 24*8=192 times  a2
				}
				{
					l3=0
					all different, one zero - comes 24*8=192 times  a3
				}
				{
					l3=l2
					largest 2 are equal - comes 12*16=192 times  b1(C)
				}
			}
			{
				l2=l1
				for(l3=l2+1; l3<L/2; ++l3)
				{
					momenta in the middle are equal - comes 12*16=192 times  b1(B)
				}
				{
					l3=L/2
					momenta in the middle equal, one L/2 -comes  12*8=192 times b2(B)
				}
				{
					l3=0
					momenta in the middle equal, one zero - 24*8=192 times b3(B)
				}
				{
					l3=l2
					largest 3 momenta are equal - comes 4*16=64 times c1(B)
				}
			}
		}
		{
			l1=l0
			for(l2=l1+1; l2<L/2; ++l2)
			{
				for(l3=l2+1; l3<L/2; ++l3)
				{
					smallest 2 are equal - comes 12*16=192 times
				}
				{
					l3=l2
					2 x 2 momenta are equal - comes 6*16=192 times
				}
			}
			{
				l2=l1
				for(l3=l2+1; l3<L/2; ++l3)
				{
					smallest 3 momenta are equal - comes 4*16=64 times
				}
				{
					//test here for zero mode if neccessary
					l3=l2
					all are equal - comes 1*16=16 time
				}
			}
		}
	}

