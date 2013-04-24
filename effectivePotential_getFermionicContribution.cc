#include "effectivePotential.h"



//compute fermionic contributions to the effective potential (only the two derivatives, since the logs are not needed but expensive)
bool effectivePotential::getFermionicContribution()
{
	if(L0==L1 && L0==L2 && L0==L3)
	{
// 		return getFermionicContribution_all_L_equal();
		return getFermionicContribution_all_L_equal_useSymmetry(); 
	}
	
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
	double one_ov_twoRho=0.5/rho;
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
					std::complex< double > w( ew*one_ov_twoRho );
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
	double one_ov_twoRho=0.5/rho;
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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




//this function uses that the eigenvalues are equal for p_mu=pi/L+x and pi/L-x
//the principle of the summation is shown in summation_principle.cc
bool effectivePotential::getFermionicContribution_all_L_equal_useSymmetry()
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
	int L=L0;
	int half_L=L/2;
	double one_ov_L(1.0/L0);
	double one_ov_twoRho=0.5/rho;
	
	U_f=0.0;
	U_fp=0.0;
	U_fpp=0.0;
	if(y_t==0.0 && y_b==0.0){ return true; }
	bool evenL= (L%2==0) ? true : false;
// 	if(!evenL)
// 	{
// 		std::cerr <<"Odd improvement not installed yet" <<std::endl;
// 		return false;
// 	}
	int loopLimit= evenL? half_L : half_L+1;
// 	std::cout <<" loopLimit: " <<loopLimit <<"   half_L: " <<half_L <<std::endl;
	
	double p0,p1,p2,p3;
	double dummy_l0_p,dummy_l1_p,dummy_l2_p; //to reduce inaccuracy by adding up a lot of doubles
	double dummy_l0_pp,dummy_l1_pp,dummy_l2_pp; //to reduce inaccuracy by adding up a lot of doubles
	std::complex< double > ew, w, z_t, z_b;
	
	for(int l0=1; l0<loopLimit; ++l0)
	{
		p0=two_PI * l0 * one_ov_L;
		dummy_l0_p=0.0; dummy_l0_pp=0.0;
		for(int l1=l0+1; l1<loopLimit; ++l1)
		{
			p1=two_PI * l1 * one_ov_L;
			dummy_l1_p=0.0; dummy_l1_pp=0.0;
			for(int l2=l1+1; l2<loopLimit; ++l2)
			{
				p2=two_PI * l2 * one_ov_L;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<loopLimit; ++l3)
				{
					p3=two_PI * l3 * one_ov_L; // a1 - 384
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 384.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 384.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=p2; // b1(C) - 192
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 192.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 192.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				if(evenL)
				{
					p3=two_PI * half_L * one_ov_L; //a2 - 192
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 192.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 192.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=0; //a3 -192
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 192.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 192.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			{
				int l2=l1;
				p2=p1;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<loopLimit; ++l3)
				{
					p3=two_PI * l3 * one_ov_L; // b1(B) - 192
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 192.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 192.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=p2; // c1(B) - 64
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 64.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 64.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				if(evenL)
				{
					p3=two_PI * half_L * one_ov_L; //b2(B) - 96
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 96.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 96.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=0; // b3(B) -96
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 96.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 96.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
			if(evenL)
			{
				p2=two_PI * half_L * one_ov_L;
				{
					p3=two_PI * half_L * one_ov_L; //b5 - 48
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 48.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 48.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=0; // a4 -96
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 96.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 96.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
			{
				p2=0;
				{
					p3=0; // b6 -48
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 48.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 48.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			dummy_l0_p+=dummy_l1_p; dummy_l0_pp+=dummy_l1_pp;
		}
		{
			int l1=l0;
			p1=p0;
			dummy_l1_p=0.0; dummy_l1_pp=0.0;
			for(int l2=l1+1; l2<loopLimit; ++l2)
			{
				p2=two_PI * l2 * one_ov_L;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<loopLimit; ++l3)
				{
					p3=two_PI * l3 * one_ov_L; // b1(A) - 192
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 192.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 192.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=p2; // d1 - 96
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 96.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 96.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				if(evenL)
				{
					p3=two_PI * half_L * one_ov_L; //b2(A) - 96
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 96.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 96.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=0; //b3(A) -96
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 96.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 96.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			{
				int l2=l1;
				p2=p1;
				dummy_l2_p=0.0; dummy_l2_pp=0.0;
				for(int l3=l2+1; l3<loopLimit; ++l3)
				{
					p3=two_PI * l3 * one_ov_L; // c1(A) - 64
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 64.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 64.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=p2; // e1  16
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 16.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 16.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"   dummy_l2_p" <<dummy_l2_p <<std::endl;
				}
				if(evenL)
				{
					p3=two_PI * half_L * one_ov_L; //c2 - 32
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 32.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 32.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=0; // c3 -32
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 32.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 32.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
			if(evenL)
			{
				p2=two_PI * half_L * one_ov_L;
				{
					p3=two_PI * half_L * one_ov_L; //d2 - 24
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 24.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 24.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=0; // b4 -48
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 48.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 48.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
			{
				p2=0;
				{
					p3=0; // d3 -24
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 24.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 24.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
		}
		if(evenL)
		{
			p1=two_PI * half_L * one_ov_L;
			{
				p2=two_PI * half_L * one_ov_L;
				{
					p3=two_PI * half_L * one_ov_L; //c4 - 8
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 8.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 8.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=0; // b7 -24
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 24.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 24.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
			{
				p2=0;
				{
					p3=0; // b8 -24
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 24.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 24.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
		}
		{
			p1=0;
			{
				p2=0;
				{
					p3=0; // c5 - 8
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 8.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 8.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				dummy_l1_p+=dummy_l2_p; dummy_l1_pp+=dummy_l2_pp;
			}
			dummy_l0_p+=dummy_l1_p; dummy_l0_pp+=dummy_l1_pp;
		}
		U_fp+=dummy_l0_p; U_fpp+=dummy_l0_pp;
	}
	dummy_l2_p=0.0; dummy_l2_pp=0.0;
	if(evenL)
	{
		p0=two_PI * half_L * one_ov_L;
		{
			p1=two_PI * half_L * one_ov_L;
			{
				p2=two_PI * half_L * one_ov_L;
				{
					p3=two_PI * half_L * one_ov_L; //e2 - 1 
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += ( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= ( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
				{
					p3=0; // c7 -4
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 4.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 4.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
			{
				p2=0;
				{
					p3=0; // d4 -6
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 6.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 6.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
		}
		{
			p1=0;
			{
				p2=0;
				{
					p3=0; // c6 - 4
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += 4.0*( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= 4.0*( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
		}
	}
	{
		p0=0;
		{
			p1=0;
			{
				p2=0;
				{
					p3=0; // e3 - 1
					ew = computeAnalyticalEigenvalue( p0,p1,p2,p3 );    w = 1.0 - ew*one_ov_twoRho;    z_t = ew + vev*y_t*w;    z_b = ew + vev*y_b*w;
					dummy_l2_p  += ( 2.0*y_t*(w/z_t).real() + 2.0*y_b*(w/z_b).real() );
					dummy_l2_pp -= ( 2.0*y_t*y_t*( w*w/(z_t*z_t)).real() + 2.0*y_b*y_b*( w*w/(z_b*z_b)).real() );
				}
			}
		}
		U_fp+=dummy_l2_p; U_fpp+=dummy_l2_pp; //NOTE different dummy here
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
	double one_ov_twoRho=0.5/rho;
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
					std::complex< double > w( ew*one_ov_twoRho );
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
