#include "effectivePotential.h"



//constructors
effectivePotential::effectivePotential(int l0, int l1, int l2, int l3, bool lowMemUsage):
L0(l0), L1(l1), L2(l2), L3(l3),
cutoff_in_GeV(-1.0),y_t(-1.0),y_b(-1.0),vev(-1.0),U_f(0.0), U_fp(0.0), U_fpp(0.0),
propagatorSumWithZeroMass(-1.0), lambda(0.0), lambda_6(0.0),
actual_HPropSum(0.0), actual_HPropSumD1(0.0), actual_HPropSumD2(0.0),
actual_GPropSum(0.0), actual_GPropSumD1(0.0), actual_GPropSumD2(0.0),
actual_BosDet(0.0), actual_BosDetD1(0.0), actual_BosDetD2(0.0),
smallest_pSquared(-1.0),lowMem(lowMemUsage)
{
	setConstants();
	if (!lowMem){ get_pSquared(); }
}

effectivePotential::effectivePotential(int l0, int l1, int l2, int l3, double cutoff, double yt, double yb, bool lowMemUsage):
L0(l0), L1(l1), L2(l2), L3(l3),
cutoff_in_GeV(cutoff),y_t(yt),y_b(yb),vev(-1.0),U_f(0.0), U_fp(0.0), U_fpp(0.0),
propagatorSumWithZeroMass(-1.0), lambda(0.0), lambda_6(0.0),
actual_HPropSum(0.0), actual_HPropSumD1(0.0), actual_HPropSumD2(0.0),
actual_GPropSum(0.0), actual_GPropSumD1(0.0), actual_GPropSumD2(0.0),
actual_BosDet(0.0), actual_BosDetD1(0.0), actual_BosDetD2(0.0),
smallest_pSquared(-1.0),lowMem(lowMemUsage)
{
	setConstants();
	vev=vev_in_GeV/cutoff_in_GeV;
	if (!lowMem){ get_pSquared(); }
	getFermionicContribution(); //for lowmem only the derivatives will be used.... will be handled in the function
}

//sets the sonsts
void effectivePotential::setConstants()
{
	r=0.5;
	rho=1.0;
	vev_in_GeV=246.0;
	N_f=1;
}




//compute fermionic contributions
bool effectivePotential::getFermionicContribution()
{
	if(L0==L1 && L0==L2 && L0==L3){ return getFermionicContribution_all_L_equal(); }
	
	if( (L0==L1 && L0==L2 && 2*L0==L3) || (L0==L1 && 2*L0==L2 && L0==L3) || (2*L0==L1 && L0==L2 && L0==L3) || (L1==L2 && L1==L3 && 2*L1==L0) )
	{
		return getFermionicContribution_Lcube_times_2L(); 
	}
	
	if(vev==-1.0 || y_t==-1.0 || y_b==-1.0)
	{
		std::cerr <<"Error, vev and/or yukawa couplings not set!" <<std::endl;
		return false;
	}
	const double two_PI(atan(1) * 8.0);
	double one_ov_L0(1.0/L0);
	double one_ov_L1(1.0/L1);
	double one_ov_L2(1.0/L2);
	double one_ov_L3(1.0/L3);
	U_f=0.0;
	U_fp=0.0;
	U_fpp=0.0;
	double dummy_l0_p,dummy_l1_p,dummy_l2_p; //to reduce inaccuracy by adding up a lot of doubles
	double dummy_l0_pp,dummy_l1_pp,dummy_l2_pp; //to reduce inaccuracy by adding up a lot of doubles
	if(y_t==0.0 && y_b==0.0){ return true; }
	for(int l0=0; l0<L0; ++l0)
	{
		double p0 = two_PI * l0 * one_ov_L0 ; //2*pi*n/L
		dummy_l0_p=0.0; dummy_l0_pp=0.0;
		for(int l1=0; l1<L1; ++l1)
		{
			double p1 = two_PI * l1 * one_ov_L1 ;
			dummy_l1_p=0.0; dummy_l1_pp=0.0;
			for(int l2=0; l2<L2; ++l2)
			{
				double p2 = two_PI * l2 * one_ov_L2 ;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=0; l3<L3; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L3 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					//U_f=\sum{log(z*z^*)} (factor of -2N_f/V comes in the end)
					//U_fp=\sum{2*y*Re[w/z]}
					//U_fpp=\sum{-2*y^2*Re[w^2/z^2]}
// 					U_f   += real( log( z_t * conj(z_t) ) + log( z_b * conj(z_b) ) ); //Never used no need
					dummy_l2_p  += 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real();
					dummy_l2_pp -= 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real();
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			dummy_l0_p+=dummy_l1_p; dummy_l0_pp+=dummy_l1_pp;
		}
		U_fp+=dummy_l0_p; U_fpp+=dummy_l0_pp;
	};
	//multiply with -2*N_f/V
// 	U_f   *= (-2.0*N_f/(L0*L1*L2*L3)); Never used, no need
	U_fp  *= (-2.0*N_f);
	U_fpp *= (-2.0*N_f);
	U_fp/=static_cast< double >(L0); U_fp/=static_cast< double >(L1); U_fp/=static_cast< double >(L2); U_fp/=static_cast< double >(L3);
	U_fpp/=static_cast< double >(L0); U_fpp/=static_cast< double >(L1); U_fpp/=static_cast< double >(L2); U_fpp/=static_cast< double >(L3);
	std::cout <<"U_fp=" <<U_fp <<"   U_fpp=" <<U_fpp <<std::endl;
	return true;
}


//computes the fermionic contribution to U' and U'' (U itself is not needed)
//this version takes advantage of the fact that all L are equal
bool effectivePotential::getFermionicContribution_all_L_equal()
{
	if(!(L0==L1 && L0==L2 && L0==L3))
	{
		std::cerr <<"Error, not all extends equal in effectivePotential::getFermionicContribution_all_L_equal()" <<std::endl;
		return false;
	}
	if(vev==-1.0 || y_t==-1.0 || y_b==-1.0)
	{
		std::cerr <<"Error, vev and/or yukawa couplings not set!" <<std::endl;
		return false;
	}
	const double two_PI(atan(1) * 8.0);
	double one_ov_L(1.0/L0);
	U_f=0.0;
	U_fp=0.0;
	U_fpp=0.0;
	if(y_t==0.0 && y_b==0.0){ return true; }
	//there are in total xx cases:
	// all momenta are different {p,q,r,s} => 4!=24 possibilities
	// 2 momenta are equal       {p,q,r,r}=> 4*3=12 possibilities NOTE catch all 3 cases, namely, that r is the smallest, medium and largest momentum!!
	// 2+2momenta are equal      {p,p,q,q}=> 4*3/2=6possibilities
	// 3 momenta are equal       {p,q,q,q}=> 4 possibilities NOTE catch all 2 cases, namely, that q is the smaller and larger momentum!!
	// all momenta are equal     {p,p,p,p}=> 1 possibiity
	
	double dummy_l0_p,dummy_l1_p,dummy_l2_p; //to reduce inaccuracy by adding up a lot of doubles
	double dummy_l0_pp,dummy_l1_pp,dummy_l2_pp; //to reduce inaccuracy by adding up a lot of doubles
	for(int l0=0; l0<L0; ++l0)
	{
		//1st: all momenta different: each kombination comes 24 (4!) times
		double p0 = two_PI * l0 * one_ov_L ; //2*pi*n/L
		dummy_l0_p=0.0; dummy_l0_pp=0.0;
		for(int l1=l0+1; l1<L1; ++l1)
		{
			double p1 = two_PI * l1 * one_ov_L ;
			dummy_l1_p=0.0; dummy_l1_pp=0.0;
			for(int l2=l1+1; l2<L2; ++l2)
			{
				double p2 = two_PI * l2 * one_ov_L ;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<L3; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					//U_f=\sum{log(z*z^*)} (factor of -2N_f/V comes in the end)
					//U_fp=\sum{2*y*Re[w/z]}
					//U_fpp=\sum{-2*y^2*Re[w^2/z^2]}
// 					U_f   += real( log( z_t * conj(z_t) ) + log( z_b * conj(z_b) ) ); //Never used no need
					dummy_l2_p  += 24.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 24.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				//2nd(a): two momenta are equal (the largest one): each combination comes 12 times
				{
					//int l3=l2; //not needed
					double p3 = p2 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 12.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 12.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			//2nd(b): two momenta are equal (the medium one): each combination comes 12 times
			{
				int l2=l1;
				double p2 = p1;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<L3; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 12.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 12.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			dummy_l0_p+=dummy_l1_p; dummy_l0_pp+=dummy_l1_pp;
		}
		//2nd(c): two momenta are equal (the smallest one): each combination comes 12 times 
		{
			//NOTE here only the case where the lowest momentum comes twice is dealt with
			int l1=l0;
			double p1 = p0 ; //meaning l1=l0
			dummy_l1_p=0.0; dummy_l1_pp=0.0;
			for(int l2=l1+1; l2<L2; ++l2)
			{
				double p2 = two_PI * l2 * one_ov_L ;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<L3; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 12.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 12.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				//3rd: 2+2 momenta are equal: comes 6 times
				{
					//int l3=l2 //not needed
					double p3 = p2 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 6.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 6.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			//4th: 3 momenta are equal: 4 possibilities NOTE only the case where the smaller is equal is treated
			{
				int l2=l1;
				double p2 = p1 ;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=0; l3<l2; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 4.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 4.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				//4th: 3 momenta are equal: 4 possibilities NOTE only the case where the larger is equal is treated
				for(int l3=l2+1; l3<L3; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 4.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 4.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				//5th: all momenta are equal: 1 possibility
				{
					//int l3=l2 //not needed
					double p3 = p2 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real();
					dummy_l2_pp -= 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real();
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			dummy_l0_p+=dummy_l1_p; dummy_l0_pp+=dummy_l1_pp;
		}
		U_fp+=dummy_l0_p; U_fpp+=dummy_l0_pp;
	}
	//multiply with -2*N_f/V
// 	U_f   *= (-2.0*N_f/(L0*L1*L2*L3)); Never used, no need
	U_fp  *= (-2.0*N_f);
	U_fpp *= (-2.0*N_f);
	U_fp/=static_cast< double >(L0); U_fp/=static_cast< double >(L1); U_fp/=static_cast< double >(L2); U_fp/=static_cast< double >(L3);
	U_fpp/=static_cast< double >(L0); U_fpp/=static_cast< double >(L1); U_fpp/=static_cast< double >(L2); U_fpp/=static_cast< double >(L3);
	std::cout <<"U_fp=" <<U_fp <<"   U_fpp=" <<U_fpp <<std::endl;
	//
	return true;
}





//computes the fermionic contribution to U' and U'' (U itself is not needed)
//this version takes advantage of the fact that the lattice looks like L^3 x 2L
bool effectivePotential::getFermionicContribution_Lcube_times_2L() //specialization for L^3x2L are equal
{
	//in this function one can abuse, that as in the L^4 case a lot of 4-touples of momenta occur multiple times 
	//for the larger extend, there are momenta, that occure nowhere, but there the sub-three-tupel can still be improved
	//in the summation l0 is interpreted as the large extend
	int L(0);
	int twoL(0);
	if(L0==L1 && L0==L2 && 2*L0==L3){ L=L0; twoL=L3 ; }
	else if(L0==L1 && 2*L0==L2 && L0==L3){ L=L0; twoL=L2 ; }
	else if(2*L0==L1 && L0==L2 && L0==L3){ L=L0; twoL=L1 ; }
	else if(L1==L2 && L1==L3 && 2*L1==L0){ L=L1; twoL=L0 ; }
	else
	{
		std::cerr <<"Error, no L^3 x 2L Lattice in effectivePotential::getFermionicContribution_Lcube_times_2L()" <<std::endl;
		return false;
	}
	
	
	if(vev==-1.0 || y_t==-1.0 || y_b==-1.0)
	{
		std::cerr <<"Error, vev and/or yukawa couplings not set!" <<std::endl;
		return false;
	}
	const double two_PI(atan(1) * 8.0);
	double one_ov_L(1.0/L);
	U_f=0.0;
	U_fp=0.0;
	U_fpp=0.0;
	if(y_t==0.0 && y_b==0.0){ return true; }
	//there are in total xx cases:
	// A (if one only considers the p_o=2*pi/L0*(2n))
	// all momenta are different {p,q,r,s} => 4!=24 possibilities
	// 2 momenta are equal       {p,q,r,r}=> 4*3=12 possibilities NOTE catch all 3 cases, namely, that r is the smallest, medium and largest momentum!!
	// 2+2momenta are equal      {p,p,q,q}=> 4*3/2=6possibilities
	// 3 momenta are equal       {p,q,q,q}=> 4 possibilities NOTE catch all 2 cases, namely, that q is the smaller and larger momentum!!
	// all momenta are equal     {p,p,p,p}=> 1 possibiity
	// B (if one only considers the p_o=2*pi/L0*(2n+1)) w in the following
	// {w,p,q,r} => 6 times
	// {w,p,q,q} => 3 times (whatch out that p can be larger and smaller then q)
	// {w,p,p,p} => 1 time
	//the case A is exactly the same as for the L^4 case, so I just copy (*crossfingers*)
	double dummy_l0_p,dummy_l1_p,dummy_l2_p; //to reduce inaccuracy by adding up a lot of doubles
	double dummy_l0_pp,dummy_l1_pp,dummy_l2_pp; //to reduce inaccuracy by adding up a lot of doubles
   //int dummyCounter=0;
	for(int l0=0; l0<L; ++l0)
	{
		//1st: all momenta different: each kombination comes 24 (4!) times
		double p0 = two_PI * l0 * one_ov_L ; //2*pi*n/L
		dummy_l0_p=0.0; dummy_l0_pp=0.0;
		for(int l1=l0+1; l1<L; ++l1)
		{
			double p1 = two_PI * l1 * one_ov_L ;
			dummy_l1_p=0.0; dummy_l1_pp=0.0;
			for(int l2=l1+1; l2<L; ++l2)
			{
				double p2 = two_PI * l2 * one_ov_L ;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<L; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					//U_f=\sum{log(z*z^*)} (factor of -2N_f/V comes in the end)
					//U_fp=\sum{2*y*Re[w/z]}
					//U_fpp=\sum{-2*y^2*Re[w^2/z^2]}
// 					U_f   += real( log( z_t * conj(z_t) ) + log( z_b * conj(z_b) ) ); //Never used no need
					dummy_l2_p  += 24.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 24.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=24;
				}
				//2nd(a): two momenta are equal (the largest one): each combination comes 12 times
				{
					//int l3=l2; //not needed
					double p3 = p2 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 12.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 12.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=12;
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			//2nd(b): two momenta are equal (the medium one): each combination comes 12 times
			{
				int l2=l1;
				double p2 = p1;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<L; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 12.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 12.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=12;
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			dummy_l0_p+=dummy_l1_p; dummy_l0_pp+=dummy_l1_pp;
		}
		//2nd(c): two momenta are equal (the smallest one): each combination comes 12 times 
		{
			int l1=l0;
			double p1 = p0 ; //meaning l1=l0
			dummy_l1_p=0.0; dummy_l1_pp=0.0;
			for(int l2=l1+1; l2<L; ++l2)
			{
				double p2 = two_PI * l2 * one_ov_L ;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<L; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 12.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 12.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=12;
				}
				//3rd: 2+2 momenta are equal: comes 6 times
				{
					//int l3=l2 //not needed
					double p3 = p2 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 6.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 6.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=6;
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			//4th: 3 momenta are equal: 4 possibilities NOTE only the case where the smaller is equal is treated
			{
				int l2=l1;
				double p2 = p1 ;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=0; l3<l2; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 4.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 4.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=4;
				}
				//4th: 3 momenta are equal: 4 possibilities NOTE only the case where the larger is equal is treated
				for(int l3=l2+1; l3<L; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 4.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 4.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=4;
				}
				//5th: all momenta are equal: 1 possibility
				{
					//int l3=l2 //not needed
					double p3 = p2 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real();
					dummy_l2_pp -= 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real();
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=1;
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			dummy_l0_p+=dummy_l1_p; dummy_l0_pp+=dummy_l1_pp;
		}
		U_fp+=dummy_l0_p; U_fpp+=dummy_l0_pp;
	}
	//now part B
	for(int l0=0; l0<L; ++l0)
	{
		double p0 = two_PI * (2*l0+1) * 0.5 * one_ov_L ; //2*pi*(n+1)/twoL
		dummy_l0_p=0.0; dummy_l0_pp=0.0;
		//6th: all three momenta are different (6 possibilities)
		for(int l1=0; l1<L; ++l1)
		{
			double p1 = two_PI * l1 * one_ov_L ;
			dummy_l1_p=0.0; dummy_l1_pp=0.0;
			for(int l2=l1+1; l2<L; ++l2)
			{
				double p2 = two_PI * l2 * one_ov_L ;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<L; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 6.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 6.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=6;
				}
				//7th(a): two momenta are equal (the larger one): each combination comes 3 times
				{
					//int l3=l2; //not needed
					double p3 = p2 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 3.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 3.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=3;
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			//7th(b): two momenta are equal (the smaller one): each combination comes 3 times
			{
				int l2=l1;
				double p2 = p1;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<L; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += 3.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 3.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=3;
				}
				//8th: all 3 momenta are equal: 1 time
				{
					//int l3=l2; //not needed
					double p3 = p2 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) ); //overlap eigenvalue nu
					std::complex< double > w( ew/(2.0*rho) );
					w=1.0-w; //w = 1 - 1/(2 rho)*nu
					std::complex< double > z_t(ew + vev*y_t*w); //z_i=nu + y*vev*w
					std::complex< double > z_b(ew + vev*y_b*w);
					dummy_l2_p  += ( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= ( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"        result_p=" <<2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() <<std::endl;
					//dummyCounter+=1;
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			dummy_l0_p+=dummy_l1_p; dummy_l0_pp+=dummy_l1_pp;
		}
		U_fp+=dummy_l0_p; U_fpp+=dummy_l0_pp;
	}
	//multiply with -2*N_f/V
// 	U_f   *= (-2.0*N_f/(L0*L1*L2*L3)); Never used, no need
	U_fp  *= (-2.0*N_f);
	U_fpp *= (-2.0*N_f);
	U_fp/=static_cast< double >(L0); U_fp/=static_cast< double >(L1); U_fp/=static_cast< double >(L2); U_fp/=static_cast< double >(L3);
	U_fpp/=static_cast< double >(L0); U_fpp/=static_cast< double >(L1); U_fpp/=static_cast< double >(L2); U_fpp/=static_cast< double >(L3);
	std::cout <<"U_fp=" <<U_fp <<"   U_fpp=" <<U_fpp <<std::endl;
	//
	return true;

}








//fills the vector with all \hat{p}^2 = 4*\sum_mu(sin(p_mu)^2)
void effectivePotential::get_pSquared()
{
	//should not be started in the lowMem case...
	if(lowMem)
	{
		std::cerr <<"Error, effectivePotential::get_pSquared(), not expected to be started with lowMem flag set!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	double p0,p1,p2,p3;
	double one_ov_L0(1.0/L0);
	double one_ov_L1(1.0/L1);
	double one_ov_L2(1.0/L2);
	double one_ov_L3(1.0/L3);
	const double PI(atan(1) * 4.0);
	for(int l0=0; l0<L0; ++l0)
	{
		p0 = sin( PI * l0 * one_ov_L0 ); //sin (p_mu /2)
		for(int l1=0; l1<L1; ++l1)
		{
			p1 = sin( PI * l1 * one_ov_L1 );
			for(int l2=0; l2<L2; ++l2)
			{
				p2 = sin( PI * l2 * one_ov_L2 );
				for(int l3=0; l3<L3; ++l3)
				{
					p3 = sin( PI * l3 * one_ov_L3 );
					pSquared.push_back( 4.0 * (p0*p0 + p1*p1 + p2*p2 + p3*p3 ) );
					
// 					std::cout <<"l0=" <<l0 <<" l1=" <<l1 <<" l2=" <<l2 <<" l3=" <<l3 <<" p0=" <<p0 <<" p1=" <<p1 <<" p2=" <<p2 <<" p3=" <<p3 <<"   p^2=" <<4.0 * (p0*p0 + p1*p1 + p2*p2 + p3*p3) <<std::endl;
				}
			}
		}
	}
	double sum(0.0);
	for(std::vector< double >::const_iterator iter=++(pSquared.begin()); iter!=pSquared.end(); ++iter)
	{
		sum+=1.0/(*iter);
	}
	propagatorSumWithZeroMass=sum/(L0*L1*L2*L3);
	//determin smallest non-zero value;
	double smallestInverse=one_ov_L0;
	if(one_ov_L1 < smallestInverse){ smallestInverse = one_ov_L1; }
	if(one_ov_L2 < smallestInverse){ smallestInverse = one_ov_L2; }
	if(one_ov_L3 < smallestInverse){ smallestInverse = one_ov_L3; }
	smallest_pSquared=4.0 * sin( PI * smallestInverse )*sin( PI * smallestInverse );
}




//return analytical eigenvalue of the overlap operator
std::complex< double > effectivePotential::computeAnalyticalEigenvalue(double p0, double p1, double p2, double p3)
{
	//computes \nu^{+} from philipp's thesis (eq 3.9)
	// \nu(p) = \rho/a + \rho/a * [i\sqrt{\tilde{p}^2} + a*r *\hat{p}^2 - \rho/a]/[\srqt{\tilde{p}^2 + (a*r *\hat{p}^2 - \rho/a)^2}]
	//a is set to 1
	//\hat{p}^2 = 1/a^2 \sum_{\mu} 4 \sin{a p_\mu/2}^2
	//\tilde{p}^2 = 1/a^2 \sum_{\mu} \sin{a p_\mu}^2
	//NOTE: in philipp's thesis it's r/2 NOTE May be an ambiguity in the definition of r...
	double p_tilde_sq=sin(p0)*sin(p0) + sin(p1)*sin(p1) + sin(p2)*sin(p2) + sin(p3)*sin(p3);
	double p_hat_sq=4.0*(sin(p0*0.5)*sin(p0*0.5) + sin(p1*0.5)*sin(p1*0.5) + sin(p2*0.5)*sin(p2*0.5) + sin(p3*0.5)*sin(p3*0.5) );
	double one_ov_denom=1.0/sqrt(p_tilde_sq + (r*p_hat_sq - rho)*(r*p_hat_sq - rho) );
	
	return std::complex< double >(rho + rho*one_ov_denom*(r*p_hat_sq - rho) , rho*one_ov_denom*sqrt(p_tilde_sq) );
}




void effectivePotential::computeBosonicDeterminantAndPropagatorSums_WithDerivatives_combined(bool updateDeterminant)
{
	//computes the bosonic determinant and the propagator sums and it's 1st and 2nd derivative wrt v
	//
	// propSum: P  = 1/( \hat{p}^2+m_0^2+z)
	//          P' = -z'/( \hat{p}^2+m_0^2+z)^2
	//          P" = [-z"*( \hat{p}^2+m_0^2+z) + 2 z'^2]/( \hat{p}^2+m_0^2+z)^3
	//det:     det = log ( \hat{p}^2+m_0^2+z)
	//         det'= z'/( \hat{p}^2+m_0^2+z)
	//         det"= [z"*( \hat{p}^2+m_0^2+z) - z'^2]/( \hat{p}^2+m_0^2+z)^2  
	if(lowMem)
	{
		std::cerr <<"Error, effectivePotential::computeBosonicDeterminantAndPropagatorSums_WithDerivatives_combined(bool updateDeterminant), not expected to be started with lowMem flag set!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(m0Squared.empty()){ std::cerr <<"Error, no value for moSquared in effectivePotential::computeBosonicDeterminantAndPropagatorSums()" <<std::endl; exit(EXIT_FAILURE); }
	
	double one_ov_V=1.0/(L0*L1*L2*L3);
	double last_m0Squared = *(--m0Squared.end());
	double z_h   = 12.0*lambda*vev*vev + 30.0*lambda_6*vev*vev*vev*vev;
	double z_hp  = 24.0*lambda*vev     + 120.0*lambda_6*vev*vev*vev;
	double z_hpp = 24.0*lambda         + 360.0*lambda_6*vev*vev;
	double z_g   = 4.0*lambda*vev*vev  + 6.0*lambda_6*vev*vev*vev*vev;
	double z_gp  = 8.0*lambda*vev      + 12.0*lambda_6*vev*vev*vev;
	double z_gpp = 8.0*lambda          + 36.0*lambda_6*vev*vev;
	
	actual_HPropSum=0.0; actual_HPropSumD1=0.0; actual_HPropSumD2=0.0;
	actual_GPropSum=0.0; actual_GPropSumD1=0.0; actual_GPropSumD2=0.0;
	actual_BosDet=0.0; actual_BosDetD1=0.0; actual_BosDetD2=0.0;
	
	if( (smallest_pSquared + last_m0Squared +z_h) <= 0.0 || (smallest_pSquared + last_m0Squared +z_g) <=0.0 )
	{
		std::cerr <<"Error, negative contribution in Propagatorsum / bosonic determinant in:" <<std::endl;
		std::cerr <<"effectivePotential::computeBosonicDeterminantAndPropagatorSums_WithDerivatives_combined(...)" <<std::endl;
		std::cerr <<"smallest_pSquared=" << smallest_pSquared <<std::endl;
		std::cerr <<"smallest_pSquared + last_m0Squared +z_h = " <<smallest_pSquared + last_m0Squared +z_h <<std::endl;
		std::cerr <<"smallest_pSquared + last_m0Squared +z_g = " <<smallest_pSquared + last_m0Squared +z_g <<std::endl;
		std::cerr <<"actual m0Squared=" << last_m0Squared <<std::endl;
		exit(EXIT_FAILURE);
	}
	
	for(std::vector< double >::const_iterator iter=++(pSquared.begin()); iter!=pSquared.end(); ++iter)
	{
		//higgs
		double propDenominator = *iter + last_m0Squared +z_h;
		actual_HPropSum   += 1.0/propDenominator;
		actual_HPropSumD1 -= z_hp/(propDenominator*propDenominator);
		actual_HPropSumD2 += (2.0*z_hp*z_hp - z_hpp*propDenominator)/(propDenominator*propDenominator*propDenominator);
		actual_BosDet     += (updateDeterminant)?0.5*log( propDenominator ):0.0;
		actual_BosDetD1   += 0.5*z_hp/propDenominator;
		actual_BosDetD2   += 0.5*(z_hpp*propDenominator - z_hp*z_hp)/(propDenominator*propDenominator);
		//goldstones
		propDenominator = *iter + last_m0Squared +z_g;
		actual_GPropSum   += 1.0/propDenominator;
		actual_GPropSumD1 -= z_gp/(propDenominator*propDenominator);
		actual_GPropSumD2 += (2.0*z_gp*z_gp - z_gpp*propDenominator)/(propDenominator*propDenominator*propDenominator);
		actual_BosDet     += (updateDeterminant)?1.5*log( propDenominator ):0.0;
		actual_BosDetD1   += 1.5*z_gp/propDenominator;
		actual_BosDetD2   += 1.5*(z_gpp*propDenominator - z_gp*z_gp)/(propDenominator*propDenominator);
	}
	actual_HPropSum   *= one_ov_V;
	actual_HPropSumD1 *= one_ov_V;
	actual_HPropSumD2 *= one_ov_V;
	actual_GPropSum   *= one_ov_V;
	actual_GPropSumD1 *= one_ov_V;
	actual_GPropSumD2 *= one_ov_V;
	actual_BosDet     *= one_ov_V;
	actual_BosDetD1   *= one_ov_V;
	actual_BosDetD2   *= one_ov_V;
	std::cout <<"last_m0Squared=" <<last_m0Squared <<std::endl; 
	std::cout <<"actual_HPropSum=" <<actual_HPropSum <<"  actual_HPropSumD1=" <<actual_HPropSumD1  <<"  actual_HPropSumD2=" <<actual_HPropSumD2 <<std::endl;
	std::cout <<"actual_GPropSum=" <<actual_GPropSum <<"  actual_GPropSumD1=" <<actual_GPropSumD1  <<"  actual_GPropSumD2=" <<actual_GPropSumD2 <<std::endl;
	std::cout <<"actual_BosDet=" <<actual_BosDet <<"  actual_BosDetD1=" <<actual_BosDetD1 <<"  actual_BosDetD2=" <<actual_BosDetD2 <<std::endl;
}

/*
void effectivePotential::computeBosonicDeterminantAndPropagatorSums_WithDerivatives_noDeterminantUpdate()
{
	//computes the 1st and second derivative of the log of the bosonic determinant and the, propagator sums and it's 1st and 2nd derivative wrt v
	//NOTE does not update the log of the bosonic determinant!!
	//
	// propSum: P  = 1/( \hat{p}^2+m_0^2+z)
	//          P' = -z'/( \hat{p}^2+m_0^2+z)^2
	//          P" = [-z"( \hat{p}^2+m_0^2+z) + 2 z'^2]/( \hat{p}^2+m_0^2+z)^3
	//det:     det = log ( \hat{p}^2+m_0^2+z)
	//         det'= z'/( \hat{p}^2+m_0^2+z)
	//         det"= [z"( \hat{p}^2+m_0^2+z) - z'^2]/( \hat{p}^2+m_0^2+z)^2  
	if(m0Squared.empty()){ std::cerr <<"Error, no value for moSquared in effectivePotential::computeBosonicDeterminantAndPropagatorSums()" <<std::endl; exit(EXIT_FAILURE); }
	
	double one_ov_V=1.0/(L0*L1*L2*L3);
	double last_m0Squared = *(--m0Squared.end());
	double z_h   = 12.0*lambda*vev*vev + 30.0*lambda_6*vev*vev*vev*vev;
	double z_hp  = 24.0*lambda*vev     + 120.0*lambda_6*vev*vev*vev;
	double z_hpp = 24.0*lambda         + 360.0*lambda_6*vev*vev;
	double z_g   = 4.0*lambda*vev*vev  + 6.0*lambda_6*vev*vev*vev*vev;
	double z_gp  = 8.0*lambda*vev      + 12.0*lambda_6*vev*vev*vev;
	double z_gpp = 8.0*lambda          + 36.0*lambda_6*vev*vev;
	
	actual_HPropSum=0.0; actual_HPropSumD1=0.0; actual_HPropSumD2=0.0;
	actual_GPropSum=0.0; actual_GPropSumD1=0.0; actual_GPropSumD2=0.0;
	actual_BosDet=0.0; actual_BosDetD1=0.0; actual_BosDetD2=0.0;
	
	if( (smallest_pSquared + last_m0Squared +z_h) <= 0.0 || (smallest_pSquared + last_m0Squared +z_g) )
	{
		std::cerr <<"Error, negative contribution in Propagatorsum / bosonic determinant in:" <<std::endl;
		std::cerr <<"effectivePotential::computeBosonicDeterminantAndPropagatorSums()" <<std::endl;
		std::cerr <<"smallest_pSquared=" << smallest_pSquared <<std::endl;
		std::cerr <<"smallest_pSquared + last_m0Squared +z_h = " <<smallest_pSquared + last_m0Squared +z_h <<std::endl;
		std::cerr <<"smallest_pSquared + last_m0Squared +z_g = " <<smallest_pSquared + last_m0Squared +z_g <<std::endl;
		std::cerr <<"actual m0Squared=" << last_m0Squared <<std::endl;
		exit(EXIT_FAILURE);
	}
	
	for(std::vector< double >::const_iterator iter=++(pSquared.begin()); iter!=pSquared.end(); ++iter)
	{
		//higgs
		double propDenominator = *iter + last_m0Squared +z_h;
		actual_HPropSum   += 1.0/propDenominator;
		actual_HPropSumD1 -= z_hp/(propDenominator*propDenominator);
		actual_HPropSumD2 += (2.0*z_hp*z_hp - z_hpp*propDenominator)/(propDenominator*propDenominator*propDenominator);
		actual_BosDetD1   += 0.5*z_hp/propDenominator;
		actual_BosDetD2   += 0.5*(z_hpp*propDenominator - z_hp*z_hp)/(propDenominator*propDenominator);
		//goldstones
		propDenominator = *iter + last_m0Squared +z_g;
		actual_GPropSum   += 1.0/propDenominator;
		actual_GPropSumD1 -= z_gp/(propDenominator*propDenominator);
		actual_GPropSumD2 += (2.0*z_gp*z_gp - z_gpp*propDenominator)/(propDenominator*propDenominator*propDenominator);
		actual_BosDetD1   += 1.5*z_gp/propDenominator;
		actual_BosDetD2   += 1.5*(z_gpp*propDenominator - z_gp*z_gp)/(propDenominator*propDenominator);
	}
	actual_HPropSum   *= one_ov_V;
	actual_HPropSumD1 *= one_ov_V;
	actual_HPropSumD2 *= one_ov_V;
	actual_GPropSum   *= one_ov_V;
	actual_GPropSumD1 *= one_ov_V;
	actual_GPropSumD2 *= one_ov_V;
	actual_BosDetD1   *= one_ov_V;
	actual_BosDetD2   *= one_ov_V;
	std::cout <<"actual_HPropSum=" <<actual_HPropSum <<"  actual_HPropSumD1=" <<actual_HPropSumD1  <<"  actual_HPropSumD2=" <<actual_HPropSumD2;
	std::cout <<"  actual_GPropSum=" <<actual_GPropSum <<"  actual_GPropSumD1=" <<actual_GPropSumD1  <<"  actual_GPropSumD2=" <<actual_GPropSumD2;
	std::cout <<"  actual_BosDetD1=" <<actual_BosDetD1 <<"  actual_BosDetD2=" <<actual_BosDetD2 <<std::endl;
}*/

void effectivePotential::computeBosonicDeterminantAndPropagatorSums_WithDerivatives_noDeterminantUpdate()
{
	computeBosonicDeterminantAndPropagatorSums_WithDerivatives_combined(false);
}

void effectivePotential::computeBosonicDeterminantAndPropagatorSums_WithDerivatives()
{
	computeBosonicDeterminantAndPropagatorSums_WithDerivatives_combined(true);
}



// void effectivePotential::computeBosonicDeterminantAndPropagatorSums_WithMassesByHand(bool updateDeterminant)
// {
// 	if(mHSquared.empty()){ std::cerr <<"Error, no value for mHSquared in effectivePotential::computeBosonicDeterminantAndPropagatorSums_WithMassesByHand()" <<std::endl; exit(EXIT_FAILURE); }
// 	actual_HPropSum=0.0; actual_HPropSumD1=0.0; actual_HPropSumD2=0.0;
// 	actual_GPropSum=0.0; actual_GPropSumD1=0.0; actual_GPropSumD2=0.0;
// 	actual_BosDet=0.0; actual_BosDetD1=0.0; actual_BosDetD2=0.0;
// 	
// 	actual_GPropSum=propagatorSumWithZeroMass;
// 	double massSquared=
// 
// }

double effectivePotential::computePropagatorSum(double massSquared)
{
// 	if( pSquared.empty() ){ get_pSquared (); } //NOTE should never occur, since it's generated in the constructor
	if( massSquared==0.0 && propagatorSumWithZeroMass > 0.0){ return propagatorSumWithZeroMass; }
	double sum(0.0);
	if (!lowMem)
	{
		for(std::vector< double >::const_iterator iter=++(pSquared.begin()); iter!=pSquared.end(); ++iter)
		{
			sum+=1.0/(*iter + massSquared);
		}
		sum/=static_cast< double >(L0); sum/=static_cast< double >(L1); sum/=static_cast< double >(L2); sum/=static_cast< double >(L3);
	}
	else if( L0==L1 && L0==L2 && L0==L3 )
	{
		sum=computePropagatorSum_lowMem_all_L_equal(massSquared);
	}
	else
	{
		double p0,p1,p2,p3;
		double one_ov_L0(1.0/L0);
		double one_ov_L1(1.0/L1);
		double one_ov_L2(1.0/L2);
		double one_ov_L3(1.0/L3);
		const double PI(atan(1) * 4.0);
		for(int l0=0; l0<L0; ++l0)
		{
			p0 = sin( PI * l0 * one_ov_L0 ); //sin (p_mu /2), p_mu=2*pi*n/L
			for(int l1=0; l1<L1; ++l1)
			{
				p1 = sin( PI * l1 * one_ov_L1 );
				for(int l2=0; l2<L2; ++l2)
				{
					p2 = sin( PI * l2 * one_ov_L2 );
					for(int l3=(l0==0 && l1==0 && l2==0 ? 1 : 0); l3<L3; ++l3) //exclude zeromode
					{
						p3 = sin( PI * l3 * one_ov_L3 );
						sum+=1.0/(4.0 * (p0*p0 + p1*p1 + p2*p2 + p3*p3 ) + massSquared);
					}
				}
			}
		}
		sum/=static_cast< double >(L0); sum/=static_cast< double >(L1); sum/=static_cast< double >(L2); sum/=static_cast< double >(L3);
	}
	if(massSquared==0.0){ propagatorSumWithZeroMass=sum; }
	
	std::cout << "effectivePotential::computePropagatorSum(" <<massSquared <<") = " <<sum <<std::endl;
	return sum; 
}

double effectivePotential::computePropagatorSum()
{
	return computePropagatorSum(0.0);
}



double effectivePotential::computePropagatorSum_lowMem_all_L_equal(double massSquared)
{
	//uses the fact, that all L are equal 
	//there are in total xx cases:
	// all momenta are different {p,q,r,s} => 4!=24 possibilities
	// 2 momenta are equal       {p,q,r,r}=> 4*3=12 possibilities NOTE catch all 3 cases, namely, that r is the smallest, medium and largest momentum!!
	// 2+2momenta are equal      {p,p,q,q}=> 4*3/2=6possibilities
	// 3 momenta are equal       {p,q,q,q}=> 4 possibilities NOTE catch all 2 cases, namely, that q is the smaller and larger momentum!!
	// all momenta are equal     {p,p,p,p}=> 1 possibiity
	//to save the computation of the sin, do it once
	
	//will also include an more accurate summation
	if(!(L0==L1 && L0==L2 && L0==L3))
	{
		std::cerr <<"Error, not all extends equal in effectivePotential::computePropagatorSum_lowMem_all_L_equal(double massSquared)" <<std::endl;
		return -1000.0;
	}
	int L=L0;
	const double PI(atan(1) * 4.0);
	double *fourTimesLatticeMomentaComponentSquared = new double[L]; // stores: 4.0*(sin (p_mu /2))^2, p_mu=2*pi*n/L
	double p0,p1,p2,p3;
	double one_ov_L(1.0/L);
	for(int l0=0; l0<L; ++l0)
	{
		p0 = sin( PI * l0 * one_ov_L ); //sin (p_mu /2), p_mu=2*pi*n/L
		fourTimesLatticeMomentaComponentSquared[l0]=4.0*p0*p0;
	}
	double sum(0.0);
	double dummy_l0, dummy_l1, dummy_l2;
	for(int l0=0; l0<L; ++l0)
	{
		p0=fourTimesLatticeMomentaComponentSquared[l0];
		dummy_l0=0.0;
		for(int l1=l0+1; l1<L; ++l1)
		{
			p1=fourTimesLatticeMomentaComponentSquared[l1];
			dummy_l1=0.0;
			for(int l2=l1+1; l2<L; ++l2)
			{
				p2=fourTimesLatticeMomentaComponentSquared[l2];
				dummy_l2=0.0;
				//1st: all momenta different: each kombination comes 24 (4!) times
				for(int l3=l2+1; l3<L; ++l3)
				{
					p3=fourTimesLatticeMomentaComponentSquared[l3];
					dummy_l2+=24.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				//2nd(a) 2 are equal (the largest): 12 times
				{
					//int l3=l2
					p3=p2;
					dummy_l2+=12.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				dummy_l1+=dummy_l2;
			}
			//2nd(b): 2 are equal (the medium one): 12 times
			{
				int l2=l1;
				p2=p1;
				dummy_l2=0.0;
				for(int l3=l2+1; l3<L; ++l3)
				{
					p3=fourTimesLatticeMomentaComponentSquared[l3];
					dummy_l2+=12.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				//4th(b): 3 momenta are equal (the larger ones): 4 times
				{
					//int l3=l2;
					p3=p2;
					dummy_l2+=4.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				dummy_l1+=dummy_l2;
			}
			dummy_l0+=dummy_l1;
		}
		//2(c):two momenta are equal (the smallest ones): 12 times
		{
			int l1=l0;
			p1=p0;
			dummy_l1=0.0;
			for(int l2=l1+1; l2<L; ++l2)
			{
				p2=fourTimesLatticeMomentaComponentSquared[l2];
				dummy_l2=0.0;
				for(int l3=l2+1; l3<L; ++l3)
				{
					p3=fourTimesLatticeMomentaComponentSquared[l3];
					dummy_l2+=12.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				//3: 2x2 momenta are equal: 6 times
				{
					//int l3=l2;
					p3=p2;
					dummy_l2+=6.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				dummy_l1+=dummy_l2;
			}
			//4th(a) 3 momenta are equal (the smaller ones): 4 times
			{
				int l2=l1;
				p2=p1;
				dummy_l2=0.0;
				for(int l3=l2+1; l3<L; ++l3)
				{
					p3=fourTimesLatticeMomentaComponentSquared[l3];
					dummy_l2+=4.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				//5th: all momenta are equal: 1 time
				//NOTE exclude zero mode!!
				if(l0!=0)
				{
					//int l3=l2;
					p3=p2;
					dummy_l2+=1.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				dummy_l1+=dummy_l2;
			}
			dummy_l0+=dummy_l1;
		}
		sum+=dummy_l0;
	}
	sum/=static_cast< double >(L0); sum/=static_cast< double >(L1); sum/=static_cast< double >(L2); sum/=static_cast< double >(L3);
	delete[] fourTimesLatticeMomentaComponentSquared;
	return sum;
}




//setting some values by hand
void effectivePotential::set_yukawa(double yt, double yb)
{
	if(yt!=y_t || yb!=y_b)
	{
		y_t=yt;
		y_b=yb;
		if(vev != -1.0){ getFermionicContribution(); }
	}
}

void effectivePotential::set_cutoff(double cutoff)
{
	if(cutoff!=cutoff_in_GeV)
	{
		cutoff_in_GeV=cutoff;
		vev=vev_in_GeV/cutoff;
		if(y_t!=-1.0 && y_b!=-1.0){ getFermionicContribution(); }
	}
}

void effectivePotential::set_lambda_6(double l_6){ lambda_6=l_6; }

void effectivePotential::set_lambda(double l){ lambda=l; }

void effectivePotential::set_rho(double new_rho){ rho=new_rho; if(y_t!=-1.0 && y_b!=-1.0 && vev !=-1.0){ getFermionicContribution(); } }

void effectivePotential::set_r(double new_r){ r = new_r; if(y_t!=-1.0 && y_b!=-1.0 && vev !=-1.0){ getFermionicContribution(); } }

void effectivePotential::set_N_f(int new_N_f){ N_f = new_N_f; if(y_t!=-1.0 && y_b!=-1.0 && vev !=-1.0){ getFermionicContribution(); } }

void effectivePotential::set_lambdaToStabilityBound() //lambda \geq -(30/48 + 24 P_G)\lambda_6.  P_G=goldstone propagator sum (without zero mode)
{
	if(lambda_6==0.0){ lambda=0.0; }
	lambda = -(30.0 / 48.0 + 24.0*computePropagatorSum(0.0))*lambda_6;
}

void effectivePotential::set_MassSquared(double new_m0Squared, double new_mHSquared)//clears the mass vectors and sets the given values as first entries
{
	m0Squared.clear(); mHSquared.clear();
	m0Squared.push_back(new_m0Squared); mHSquared.push_back(new_mHSquared);
}

void effectivePotential::resetConstants(double new_rho, double new_r, int new_N_f)
{
	if(rho!=new_rho || r!=new_r || N_f!=new_N_f)
	{
		rho=new_rho;
		r=new_r;
		N_f=new_N_f;
		if(vev != -1.0 && y_t!=-1.0 && y_b!=-1.0){ getFermionicContribution(); }
	}
}


double effectivePotential::get_y_t(){ return y_t; }

double effectivePotential::get_y_b(){ return y_b; } 

double effectivePotential::get_cutoff_in_GeV(){ return cutoff_in_GeV; }

double effectivePotential::get_vev(){ return vev; }

double effectivePotential::get_lambda_6(){ return lambda_6; }

double effectivePotential::get_lambda(){ return lambda; }

double effectivePotential::get_U_f(){ return U_f; }

double effectivePotential::get_U_fp(){ return U_fp; }

double effectivePotential::get_U_fpp(){ return U_fpp; }

std::vector< double > effectivePotential::get_m0Squared_sequence(){ return m0Squared; }

std::vector< double > effectivePotential::get_mHSquared_sequence(){ return mHSquared; } 

double effectivePotential::get_m0Squared(){ if(!(m0Squared.empty())){ return *(m0Squared.rbegin()); } else{ return 0.0; } }

double effectivePotential::get_mHSquared(){ if(!(mHSquared.empty())){ return *(mHSquared.rbegin()); } else{ return 0.0; } }

double effectivePotential::get_mH_in_GeV(){ if(!(mHSquared.empty())){ return cutoff_in_GeV*sqrt(*(mHSquared.rbegin())); } else{ return 0.0; } }


void effectivePotential::initializeTreeLevel()
{
	if(vev==-1.0 || y_t ==-1.0 || y_b==-1.0)
	{
		std::cerr <<"Error, vev, y_t or y_b not defined in effectivePotential::initializeTreeLevel()" <<std::endl;
		exit(EXIT_FAILURE);
	}
	//clean vecotrs if not empty
	m0Squared.clear();
	mHSquared.clear();
	double m0Squared_tree=-U_fp/vev - 4.0*vev*vev*lambda -6.0*vev*vev*vev*vev*lambda_6;
// 	if(m0Squared_tree < 0.0){ m0Squared.push_back(0.0); }//NOTE why would you do this?? m0^2 can be negative
// 	else{ m0Squared.push_back(m0Squared_tree); }
	m0Squared.push_back(m0Squared_tree);	
	double mHSquared_tree=U_fpp + m0Squared_tree + 12.0*vev*vev*lambda + 30.0*vev*vev*vev*vev*lambda_6;
// 	mHSquared_tree=0.0;
	if(mHSquared_tree < 0.0){ mHSquared.push_back(0.0); } //might be reasonable at treelevel
	else{ mHSquared.push_back(mHSquared_tree); }
}




double effectivePotential::iterateMass_withM0inPropSumsAndBosonicDet()
{
	//NOTE not lowMem support ;-) does not work anyway
	if(lowMem)
	{
		std::cerr <<"Error, effectivePotential::iterateMass_withM0inPropSumsAndBosonicDet() not expected to be started with lowMem flag set!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(m0Squared.empty() || mHSquared.empty()){ std::cerr <<"no startvalues for iteration" <<std::endl; initializeTreeLevel(); }
	
	//
	computeBosonicDeterminantAndPropagatorSums_WithDerivatives_noDeterminantUpdate();
	double last_m0Squared = *(--m0Squared.end());
	// U = U_f + 0.5*m_0^2*v^2 + lambda*v^4 + lambda_6*v^6 + BosDet +
	//   + v^2*[ 6*lambda*(P_G + P_H) + lambda_6*(45*P_G^2 + 54*P_G*P_H + 45*P_H^2) ] +
	//   + v^4*[ lambda_6*(15*P_H + 9*P_G) ]
	// U'= U_f' + m_0^2*v + 4*lambda*v^3 +6*lambda_6*v^5 +BosDet' +
	//   + 2*v*[ 6*lambda*(P_G + P_H) + lambda_6*(45*P_G^2 + 54*P_G*P_H + 45*P_H^2) ] + v^2*[ 6*lambda*(P_G' + P_H') + lambda_6*( 90*P_G*P_G' + 54*(P_G'*P_H + P_G*P_H') + 90*P_H*P_H') ]
	//   + 4*v^3*[ lambda_6*(15*P_H + 9*P_G) ] + v^4*[ lambda_6*(15*P_H' + 9*P_G') ]
	// U"= U_F" + m_0^2 + 12*lambda*v^2 + 30*lambda_6*v^4 +BosDet" +
	//   + 2*[ 6*lambda*(P_G + P_H) + lambda_6*(45*P_G^2 + 54*P_G*P_H + 45*P_H^2) ] +2*v*[ 6*lambda*(P_G' + P_H') + lambda_6*( 90*P_G*P_G' + 54*(P_G'*P_H + P_G*P_H') + 90*P_H*P_H') ] +
	//   + 2*v*[ 6*lambda*(P_G' + P_H') + lambda_6*( 90*P_G*P_G' + 54*(P_G'*P_H + P_G*P_H') + 90*P_H*P_H') ] + 
	//   + v^2*[ 6*lambda*(P_G" + P_H") + lambda_6*( 90*(P_G*P_G" + P_G'^2) + 54*(P_G"*P_H + 2*P_G'*P_H' + P_G*P_H") + 90*(P_H*P_H" + P_H'^2) ) ]
	//   + 12*v^2*[ lambda_6*(15*P_H + 9*P_G) ] + 4*v^3*[ lambda_6*(15*P_H' + 9*P_G') ]
	//   + 4*v^3*[ lambda_6*(15*P_H' + 9*P_G') ] + v^4*[ lambda_6*(15*P_H" + 9*P_G") ]
	//
	//in short:
	//U = U_tree + 0.5*m_0^2*v^2 + v^2*A + v^4*B
	//U'= U_tree' + m_0^2*v + 2*v*A + v^2*A' + 4*v^3*B + v^4*B'  =>  m_0^2=-(U_tree'/v +2A + v*A' +4*v^2*B + v^3*B')
	//U"= U_tree" + m_0^2 + 2*A + 2*v*A' +2*v*A' + v^2*A" + 12*v^2*B + 4*v^3*B' + 4*v^3*B' + v^4*B"
	//  = U_tree" + m_0^2 + 2*A + 4*v*A' + v^2*A" + 12*v^2*B + 8*v^3*B' + v^4*B"
	double U_tree_p = U_fp + 4.0*lambda*vev*vev*vev + 6.0*lambda_6*vev*vev*vev*vev*vev + actual_BosDetD1;
	double U_tree_pp = U_fpp + 12.0*lambda*vev*vev + 30.0*lambda_6*vev*vev*vev*vev + actual_BosDetD2;
	double A = 6.0*lambda*( actual_HPropSum + actual_GPropSum ) 
	         + lambda_6*( 45.0*actual_GPropSum*actual_GPropSum 
	                    + 54.0*actual_GPropSum*actual_HPropSum 
	                    + 45.0*actual_HPropSum*actual_HPropSum );
	double A_p = 6.0*lambda*( actual_HPropSumD1 + actual_GPropSumD1 )
	           + lambda_6*( 90.0*actual_GPropSum*actual_GPropSumD1 
	                      + 54.0*( actual_GPropSumD1*actual_HPropSum + actual_GPropSum*actual_HPropSumD1 )
	                      + 90.0*actual_HPropSum*actual_HPropSumD1 );
	double A_pp = 6.0*lambda*( actual_HPropSumD2 + actual_GPropSumD2 )
	            + lambda_6*( 90.0*( actual_GPropSum*actual_GPropSumD2 + actual_GPropSumD1*actual_GPropSumD1 )
	                       + 54.0*( actual_GPropSumD2*actual_HPropSum + 2.0*actual_GPropSumD1*actual_HPropSumD1 + actual_GPropSum*actual_HPropSumD2)
	                       + 90.0*( actual_HPropSum*actual_HPropSumD2 + actual_HPropSumD1*actual_HPropSumD1 ) );
	double B = lambda_6*( 15.0*actual_HPropSum + 9.0*actual_GPropSum );
	double B_p = lambda_6*( 15.0*actual_HPropSumD1 + 9.0*actual_GPropSumD1 );
	double B_pp = lambda_6*( 15.0*actual_HPropSumD2 + 9.0*actual_GPropSumD2 );
	
	double new_m0Squared = -(U_tree_p/vev + 2.0*A + vev*A_p + 4.0*vev*vev*B + vev*vev*vev*B_p);
	double new_mHSquared = U_tree_pp + new_m0Squared + 2.0*A + 4.0*vev*A_p + vev*vev*A_pp + 12.0*vev*vev*B + 8.0*vev*vev*vev*B_p + vev*vev*vev*vev*B_pp;
	m0Squared.push_back(new_m0Squared);
	mHSquared.push_back(new_mHSquared);
	std::cout <<"lambda=" <<lambda <<"  lambda_6=" <<lambda_6 <<"  U_fp" <<U_fp <<"  U_fpp" <<U_fpp <<std::endl;
	std::cout <<"U_tree_p=" <<U_tree_p <<"  U_tree_pp=" <<U_tree_pp;
	std::cout <<"  A=" <<A <<"  A_p=" <<A_p <<"  A_pp=" <<A_pp;
	std::cout <<"  B=" <<B <<"  B_p=" <<B_p <<"  B_pp=" <<B_pp <<"  mHSquared=" <<new_mHSquared <<std::endl; 
	return( std::abs(new_m0Squared - last_m0Squared) );
}





double effectivePotential::iterateMassDetermination_withMassInPropSumByHand()
{
	if(m0Squared.empty() || mHSquared.empty()){ std::cerr <<"no startvalues for iteration" <<std::endl; initializeTreeLevel(); }
	double last_mHSquared = *(--mHSquared.end());
	double last_m0Squared = *(--m0Squared.end());
// 	std::cout <<"last_mHSquared:  " <<last_mHSquared  <<std::endl;
	double GoldstonePropagatorSum=computePropagatorSum(0.0);
	double HiggsPropagatorSum=computePropagatorSum(last_mHSquared);
	//new m0
	double new_m0Squared=-U_fp/vev -4.0*vev*vev*lambda -6.0*vev*vev*vev*vev*lambda_6; //treelevel
// 	std::cout <<"new_m0Squared 1 " <<new_m0Squared <<std::endl;
	new_m0Squared+= -12.0*lambda*(HiggsPropagatorSum + GoldstonePropagatorSum);
// 	std::cout <<"new_m0Squared 2 " <<new_m0Squared <<std::endl;
	new_m0Squared+= -lambda_6*(90.0*HiggsPropagatorSum*HiggsPropagatorSum + 108.0*GoldstonePropagatorSum*HiggsPropagatorSum + 90.0*GoldstonePropagatorSum*GoldstonePropagatorSum);
// 	std::cout <<"new_m0Squared 3 " <<new_m0Squared <<std::endl;
	new_m0Squared+= -vev*vev*lambda_6*(60.0*HiggsPropagatorSum + 36.0*GoldstonePropagatorSum);
// 	std::cout <<"new_m0Squared 4 " <<new_m0Squared <<std::endl;
	double new_mHSquared=U_fpp + new_m0Squared + 12.0*vev*vev*lambda + 30.0*vev*vev*vev*vev*lambda_6; //treelevel
	new_mHSquared += 12.0*lambda*(HiggsPropagatorSum + GoldstonePropagatorSum);
	new_mHSquared += lambda_6*(90.0*HiggsPropagatorSum*HiggsPropagatorSum + 108.0*GoldstonePropagatorSum*HiggsPropagatorSum + 90.0*GoldstonePropagatorSum*GoldstonePropagatorSum);
	new_mHSquared += vev*vev*lambda_6*(180.0*HiggsPropagatorSum + 108.0*GoldstonePropagatorSum);
	m0Squared.push_back(new_m0Squared);
	mHSquared.push_back(new_mHSquared);
	return (std::abs(new_m0Squared - last_m0Squared) + std::abs(new_mHSquared - last_mHSquared) );
}




bool effectivePotential::fillVectorWithEigenvaluesAndPrefac( std::vector< std::pair< std::complex< double >, std::complex< double > > > &toBeFilled )
{
	//not with lowMem feature!!
	if(lowMem)
	{
		std::cerr <<"Error, effectivePotential::fillVectorWithEigenvaluesAndPrefac, not expected to be started with lowMem flag set!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(y_t==-1.0 || y_b==-1.0)
	{
		std::cerr <<"Error not set!" <<std::endl;
		return false;
	}
	const double two_PI(atan(1) * 8.0);
	double one_ov_L0(1.0/L0);
	double one_ov_L1(1.0/L1);
	double one_ov_L2(1.0/L2);
	double one_ov_L3(1.0/L3);
	for(int l0=0; l0<L0; ++l0)
	{
		double p0 = two_PI * l0 * one_ov_L0 ; //2*pi*n/L
		for(int l1=0; l1<L1; ++l1)
		{
			double p1 = two_PI * l1 * one_ov_L1 ;
			for(int l2=0; l2<L2; ++l2)
			{
				double p2 = two_PI * l2 * one_ov_L2 ;
				for(int l3=0; l3<L3; ++l3)
				{
					double p3 = two_PI * l3 * one_ov_L3 ;
					std::complex< double > ew( computeAnalyticalEigenvalue( p0,p1,p2,p3 ) );
					std::complex< double > w(1.0-ew/(2.0*rho) );
					toBeFilled.push_back(std::make_pair(ew, w) );
				}
			}
		}
	}
	return true;
}





bool effectivePotential::plotPotentialToStream_withMassInPropSumByHand(double min, double max, double step, std::ostream &output, bool withInfo)
{
	//U=U_f + phi^2*[m0^2/2 + 6*lambda*(P_H+P_G) + lambda_6*(45*P_H^2+54*P_H*P_G+45*P_G^2)] + phi^4[lambda + lambda_6*(15*P_h+9*P_G)] + phi^6*lambda_6
	//U=U_f + phi^2*A + phi^4*B + phi^6*lambda_6
	if(lowMem)
	{
		std::cerr <<"Error, effectivePotential::scanPotential_firstDerivative_withMassInPropSumByHand, not expected to be started with lowMem flag set!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(m0Squared.empty() || mHSquared.empty()){ std::cerr <<"no masses available" <<std::endl; return false; }
	if((max-min)/step > maxCount || (max-min) < 0.0 || step <=0.0  )
	{
		std::cerr <<"Error, inconstistent range or to many steps. maxCount in effectivePotential::plotPotentialToStream is " <<maxCount <<std::endl;
		return false;
	}
	//increase precision, but store the old one
	std::streamsize oldPrec=output.precision();
	output.precision(15);
	double last_m0Squared=*(m0Squared.rbegin());
	double last_mHSquared=*(mHSquared.rbegin());
	double HiggsPropagatorSum=computePropagatorSum(last_mHSquared);
	std::vector< std::pair< std::complex< double >, std::complex< double > > > eigenvaluesAndPrefactor;// contains: < nu, (1-1/(2*rho))nu >, nu are the eigenvalues of the overlap operator
	bool withFermions(y_t!=0.0 || y_b!=0.0);
	if( withFermions)
	{
		if(!fillVectorWithEigenvaluesAndPrefac(eigenvaluesAndPrefactor))
		{
			std::cerr <<"Error, generating Eigenvalues" <<std::endl;
			return false;
		}
	}
	double fermionicPrefactor(-2.0*N_f/(L0*L1*L2*L3));
	//A=[m0^2/2 + 6*lambda*(P_H+P_G) + lambda_6*(45*P_H^2+54*P_H*P_G+45*P_G^2)]
	//B=[lambda + lambda_6*(15*P_h+9*P_G)]
	double A( 0.5*last_m0Squared + 6.0*lambda*(HiggsPropagatorSum + propagatorSumWithZeroMass)
	      + lambda_6*(45.0*HiggsPropagatorSum*HiggsPropagatorSum + 54.0*HiggsPropagatorSum*propagatorSumWithZeroMass + 45.0*propagatorSumWithZeroMass*propagatorSumWithZeroMass) );
	double B( lambda + lambda_6*(15.0*HiggsPropagatorSum + 9.0*propagatorSumWithZeroMass) );
	double phi(min);
	double U(0.0), U_p(0.0), U_pp(0.0);
	if(withInfo)
	{
		output <<"# Values for the effective potential to first order in lambda and lambda_6" <<std::endl;
		output <<"# values used for the potential are:" <<std::endl;
		output <<"# Cutoff in GeV: " <<cutoff_in_GeV  <<std::endl;
		output <<"# vev: " <<vev  <<std::endl;
		output <<"# y_t: " <<y_t  <<std::endl;
		output <<"# y_b: " <<y_b  <<std::endl;
		output <<"# lambda_6: " <<lambda_6 <<std::endl;
		output <<"# lambda: " <<lambda <<std::endl;
		output <<"# m_0^2: " <<last_m0Squared <<std::endl;
		output <<"# m_H^2: " <<last_mHSquared <<std::endl;
		output <<"# m_H in GeV: " <<sqrt(last_mHSquared)*cutoff_in_GeV <<std::endl;
		output <<"# format is: phi   U[phi]   U'[phi]   U''[phi]" <<std::endl;
	}
	while(phi <= max)
	{
		U=0.0;
		U_p=0.0;
		U_pp=0.0;
		//fermionic contribution
		if(withFermions)
		{
			for(std::vector< std::pair< std::complex< double >, std::complex< double > > >::const_iterator iter=eigenvaluesAndPrefactor.begin(); iter!=eigenvaluesAndPrefactor.end(); ++iter)
			{
				//a clearer way can be found in effectivePotential::getFermionicContribution()
				std::complex< double > z_t(iter->first + y_t*phi*iter->second);
				std::complex< double > z_b(iter->first + y_b*phi*iter->second);
				U += log( real( z_t * conj(z_t) ) ) + log( real( z_b * conj(z_b) ) );
				U_p += 2.0*y_t*real( iter->second / z_t ) + 2.0*y_b*real( iter->second / z_b );
				U_pp -= 2.0*y_t*y_t*real( iter->second * iter->second /(z_t * z_t) ) + 2.0*y_b*y_b*real( iter->second * iter->second /(z_b * z_b) );
			}
			U*=fermionicPrefactor;
			U_p*=fermionicPrefactor;
			U_pp*=fermionicPrefactor;
		}
		U    +=phi*phi*A + phi*phi*phi*phi*B + phi*phi*phi*phi*phi*phi*lambda_6;
		U_p  +=2.0*phi*A + 4.0*phi*phi*phi*B + 6.0*phi*phi*phi*phi*phi*lambda_6;
		U_pp +=2.0*A     + 12.0*phi*phi*B    + 30.0*phi*phi*phi*phi*lambda_6;
		//increase precision, but store the old one
		if(!(output <<phi <<" " <<U <<" " <<U_p <<" " <<U_pp <<std::endl))
		{
			std::cerr <<"Error during output in effectivePotential::plotPotentialToStream" <<std::endl;
			return false;
		}
		phi+=step;
	}
	output.precision(oldPrec);
	return true;
}



bool effectivePotential::scanPotential_firstDerivative_withMassInPropSumByHand(double min, double max, double step, std::map< double, double > &result_U_p)
{
	// computes the first derivative of the potential in a range of v from min to max in steps of ...(guess!)
	//U_p = U_fp[v] + m0Squared*v + 4*lambda*v^3 + 6*lambda_6*v^5 + 
	//    + 2*v*(6*lambda*(P_H + P_G) + lambda_6*(45*P_G^2 + 54*P_G*P_H + 45*P_H^2) ) + 
	//    + 4*v^3*lambda_6*(15*P_H + 9*P_G)
	//    == U_fp[v] + v*A + v^3*B + 6*lambda_6*v^5
	//
	//U_fp[v] = -2*N_f/V * 2*y_t * sum_p[ re((1-0.5*rho*nu)/(nu+y_t*v*(1-0.5*rho*nu))) ] + (same with y_b) 
	
	if(lowMem)
	{
		std::cerr <<"Error, effectivePotential::scanPotential_firstDerivative_withMassInPropSumByHand, not expected to be started with lowMem flag set!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(m0Squared.empty() || mHSquared.empty()){ std::cerr <<"no masses available" <<std::endl; return false; }
	
	result_U_p.clear();
	
	//consistency
	if((max-min)/step > maxCount || (max-min) < 0.0 || step <=0.0  )
	{
		std::cerr <<"Error, inconstistent range or to many steps. maxCount in effectivePotential::plotPotentialToStream is " <<maxCount <<std::endl;
		return false;
	}
	//
	double last_m0Squared=*(m0Squared.rbegin());
	double last_mHSquared=*(mHSquared.rbegin());
	double HiggsPropagatorSum=computePropagatorSum(last_mHSquared);
	
	// contains: < nu, (1-1/(2*rho))nu >, nu are the eigenvalues of the overlap operator
	std::vector< std::pair< std::complex< double >, std::complex< double > > > eigenvaluesAndPrefactor;
	bool withFermions(y_t!=0.0 || y_b!=0.0);
	if( withFermions)
	{
		if(!fillVectorWithEigenvaluesAndPrefac(eigenvaluesAndPrefactor))
		{
			std::cerr <<"Error, generating Eigenvalues" <<std::endl;
			return false;
		}
	}
	double fermionicPrefactor(-2.0*N_f/(L0*L1*L2*L3));
	//A=2*(6*lambda*(P_H + P_G) + lambda_6*(45*P_G^2 + 54*P_G*P_H + 45*P_H^2) ) + m0Squared
	//B=4*lambda_6*(15*P_H + 9*P_G) + 4*lambda
	double A( last_m0Squared + 12.0*lambda*(HiggsPropagatorSum + propagatorSumWithZeroMass)
	      + lambda_6*(90.0*HiggsPropagatorSum*HiggsPropagatorSum + 108.0*HiggsPropagatorSum*propagatorSumWithZeroMass + 90.0*propagatorSumWithZeroMass*propagatorSumWithZeroMass) );
	double B(4.0*lambda + lambda_6*(60.0*HiggsPropagatorSum + 36.0*propagatorSumWithZeroMass) );
	double phi(min);
	double U_p(0.0);
	std::map< double, double >::iterator hintIter=result_U_p.begin(); //faster insertion
	while(phi <= max)
	{
		U_p=0.0;
		//fermionic contribution
		if(withFermions)
		{
			for(std::vector< std::pair< std::complex< double >, std::complex< double > > >::const_iterator iter=eigenvaluesAndPrefactor.begin(); iter!=eigenvaluesAndPrefactor.end(); ++iter)
			{
				//a clearer way can be found in effectivePotential::getFermionicContribution()
				std::complex< double > z_t(iter->first + y_t*phi*iter->second);
				std::complex< double > z_b(iter->first + y_b*phi*iter->second);
				U_p += 2.0*y_t*real( iter->second / z_t ) + 2.0*y_b*real( iter->second / z_b );
			}
			U_p*=fermionicPrefactor;
		}
		U_p  += phi*A + phi*phi*phi*B + 6.0*phi*phi*phi*phi*phi*lambda_6;
		hintIter=result_U_p.insert(hintIter, std::make_pair( phi, U_p) ); 
		phi+=step;
	}
	return true;

}


bool effectivePotential::scanPotential_secondDerivative_withMassInPropSumByHand(double min, double max, double step, std::map< double, double > &result_U_pp)
{
	// computes the second derivative of the potential in a range of v from min to max in steps of ...(guess!)
	//U_pp = U_fpp[v] + m0Squared + 12*lambda*v^2 + 30*lambda_6*v^4 + 
	//    + 2*(6*lambda*(P_H + P_G) + lambda_6*(45*P_G^2 + 54*P_G*P_H + 45*P_H^2) ) + 
	//    + 3*4*v^2*lambda_6*(15*P_H + 9*P_G)
	//    == U_fpp[v] + A + 3*v^2*B + 30*lambda_6*v^4
	//
	//    A=2*(6*lambda*(P_H + P_G) + lambda_6*(45*P_G^2 + 54*P_G*P_H + 45*P_H^2) ) + m0Squared
	//    B=4*lambda_6*(15*P_H + 9*P_G) + 12*lambda
	//
	//U_fpp[v] = -2*N_f/V * (-2*y_t^2 * sum_p[ re( (1-0.5*rho*nu)^2/(nu+y_t*v*(1-0.5*rho*nu))^2 ) ] + (same with y_b) 
	if(lowMem)
	{
		std::cerr <<"Error, effectivePotential::scanPotential_secondDerivative_withMassInPropSumByHand, not expected to be started with lowMem flag set!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(m0Squared.empty() || mHSquared.empty()){ std::cerr <<"no masses available" <<std::endl; return false; }
	
	result_U_pp.clear();
	
	//consistency
	if((max-min)/step > maxCount || (max-min) < 0.0 || step <=0.0  )
	{
		std::cerr <<"Error, inconstistent range or to many steps. maxCount in effectivePotential::plotPotentialToStream is " <<maxCount <<std::endl;
		return false;
	}
	double last_m0Squared=*(m0Squared.rbegin());
	double last_mHSquared=*(mHSquared.rbegin());
	double HiggsPropagatorSum=computePropagatorSum(last_mHSquared);
	
	// contains: < nu, (1-1/(2*rho))nu >, nu are the eigenvalues of the overlap operator
	std::vector< std::pair< std::complex< double >, std::complex< double > > > eigenvaluesAndPrefactor;
	bool withFermions(y_t!=0.0 || y_b!=0.0);
	if( withFermions)
	{
		if(!fillVectorWithEigenvaluesAndPrefac(eigenvaluesAndPrefactor))
		{
			std::cerr <<"Error, generating Eigenvalues" <<std::endl;
			return false;
		}
	}
	double fermionicPrefactor(-2.0*N_f/(L0*L1*L2*L3));
	//A=2*(6*lambda*(P_H + P_G) + lambda_6*(45*P_G^2 + 54*P_G*P_H + 45*P_H^2) ) + m0Squared
	//B=4*lambda_6*(15*P_H + 9*P_G) + 4*lambda
	double A( last_m0Squared + 12.0*lambda*(HiggsPropagatorSum + propagatorSumWithZeroMass)
	      + lambda_6*(90.0*HiggsPropagatorSum*HiggsPropagatorSum + 108.0*HiggsPropagatorSum*propagatorSumWithZeroMass + 90.0*propagatorSumWithZeroMass*propagatorSumWithZeroMass) );
	double B(4.0*lambda + lambda_6*(60.0*HiggsPropagatorSum + 36.0*propagatorSumWithZeroMass) );
	double phi(min);
	double U_pp(0.0);
	std::map< double, double >::iterator hintIter=result_U_pp.begin(); //faster insertion
	while(phi <= max)
	{
		U_pp=0.0;
		//fermionic contribution
		if(withFermions)
		{
			for(std::vector< std::pair< std::complex< double >, std::complex< double > > >::const_iterator iter=eigenvaluesAndPrefactor.begin(); iter!=eigenvaluesAndPrefactor.end(); ++iter)
			{
				//a clearer way can be found in effectivePotential::getFermionicContribution()
				std::complex< double > z_t(iter->first + y_t*phi*iter->second);
				std::complex< double > z_b(iter->first + y_b*phi*iter->second);
				U_pp -= 2.0*y_t*y_t*real( (iter->second * iter->second)/(z_t*z_t) ) + 2.0*y_b*y_b*real( (iter->second * iter->second)/(z_b*z_b) );
			}
			U_pp*=fermionicPrefactor;
		}
		U_pp  += A + 3.0*phi*phi*B + 30.0*phi*phi*phi*phi*lambda_6;
		hintIter=result_U_pp.insert(hintIter, std::make_pair( phi, U_pp) ); 
		phi+=step;
	}
	return true;

}	
	
