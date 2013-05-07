#include "effectivePotential.h"


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
// 		sum=computePropagatorSum_lowMem_all_L_equal(massSquared);
		sum=computePropagatorSum_lowMem_all_L_equal_useSymmetry(massSquared);
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
// 						std::cout <<"l0=" <<l0 <<"  l1=" <<l1 <<"  l2=" <<l2 <<"  l3=" <<l3
// 						     <<"   p0=" <<p0 <<"  p1=" <<p1 <<"  p2=" <<p2 <<"  p3=" <<p3
// 						     <<"1.0/(4.0 * (p0*p0 + p1*p1 + p2*p2 + p3*p3 ) + " <<massSquared <<")="
// 						     <<1.0/(4.0 * (p0*p0 + p1*p1 + p2*p2 + p3*p3 ) + massSquared) <<std::endl;
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


//this function uses that the contributions are equal for p_mu=(pi+x)/L and (pi-x)/L
//the principle of the summation is shown in summation_principle.cc
double effectivePotential::computePropagatorSum_lowMem_all_L_equal_useSymmetry(double massSquared)
{
	if(!(L0==L1 && L0==L2 && L0==L3))
	{
		std::cerr <<"Error, not all extends equal in effectivePotential::computePropagatorSum_lowMem_all_L_equal(double massSquared)" <<std::endl;
		return -1000.0;
	}
	int L=L0;
	const double PI(atan(1) * 4.0);
	int half_L=L/2;
	bool evenL= (L%2==0) ? true : false;
	int loopLimit= evenL? half_L : half_L+1;
	
	
	double *fourTimesLatticeMomentaComponentSquared = new double[half_L]; // stores: 4.0*(sin (p_mu /2))^2, p_mu=2*pi*n/L
	double p0,p1,p2,p3;
	double one_ov_L(1.0/L);
	for(int l0=0; l0<=half_L; ++l0)
	{
		p0 = sin( PI * l0 * one_ov_L ); //sin (p_mu /2), p_mu=2*pi*n/L
		fourTimesLatticeMomentaComponentSquared[l0]=4.0*p0*p0;
	}
	double sum(0.0);
	double dummy_l0, dummy_l1, dummy_l2;
	
	for(int l0=1; l0<loopLimit; ++l0)
	{
		p0=fourTimesLatticeMomentaComponentSquared[l0];
		dummy_l0=0.0;
		for(int l1=l0+1; l1<loopLimit; ++l1)
		{
			p1=fourTimesLatticeMomentaComponentSquared[l1];
			dummy_l1=0.0;
			for(int l2=l1+1; l2<loopLimit; ++l2)
			{
				p2=fourTimesLatticeMomentaComponentSquared[l2];
				dummy_l2=0.0;
				for(int l3=l2+1; l3<loopLimit; ++l3)
				{
					p3=fourTimesLatticeMomentaComponentSquared[l3]; // a1 - 384
					dummy_l2 += 384.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=p2; // b1(C) - 192
					dummy_l2 += 192.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				if(evenL)
				{
					p3=fourTimesLatticeMomentaComponentSquared[half_L]; //a2 - 192
					dummy_l2 += 192.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=0; //a3 -192
					dummy_l2 += 192.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				dummy_l1+=dummy_l2;
			}
			{
				int l2=l1;
				p2=p1;
				dummy_l2=0.0;
				for(int l3=l2+1; l3<loopLimit; ++l3)
				{
					p3=fourTimesLatticeMomentaComponentSquared[l3]; // b1(B) - 192
					dummy_l2 += 192.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=p2; // c1(B) - 64
					dummy_l2 += 64.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				if(evenL)
				{
					p3=fourTimesLatticeMomentaComponentSquared[half_L]; //b2(B) - 96
					dummy_l2 += 96.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=0; // b3(B) -96
					dummy_l2 += 96.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
			if(evenL)
			{
				p2=fourTimesLatticeMomentaComponentSquared[half_L];
				{
					p3=fourTimesLatticeMomentaComponentSquared[half_L]; //b5 - 48
					dummy_l2 += 48.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=0; // a4 -96
					dummy_l2 += 96.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
			{
				p2=0;
				{
					p3=0; // b6 -48
					dummy_l2 += 48.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				dummy_l1+=dummy_l2;
			}
			dummy_l0+=dummy_l1;
		}
		{
			int l1=l0;
			p1=p0;
			dummy_l1=0.0;
			for(int l2=l1+1; l2<loopLimit; ++l2)
			{
				p2=fourTimesLatticeMomentaComponentSquared[l2];
				dummy_l2=0.0;
				for(int l3=l2+1; l3<loopLimit; ++l3)
				{
					p3=fourTimesLatticeMomentaComponentSquared[l3]; // b1(A) - 192
					dummy_l2 += 192.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=p2; // d1 - 96
					dummy_l2 += 96.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				if(evenL)
				{
					p3=fourTimesLatticeMomentaComponentSquared[half_L]; //b2(A) - 96
					dummy_l2 += 96.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=0; //b3(A) -96
					dummy_l2 += 96.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				dummy_l1+=dummy_l2;
			}
			{
				int l2=l1;
				p2=p1;
				dummy_l2=0.0;
				for(int l3=l2+1; l3<loopLimit; ++l3)
				{
					p3=fourTimesLatticeMomentaComponentSquared[l3]; // c1(A) - 64
					dummy_l2 += 64.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=p2; // e1  16
					dummy_l2 += 16.0/( p0 + p1 + p2 + p3 + massSquared);
					//std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<"   dummy_l2_p" <<dummy_l2_p <<std::endl;
				}
				if(evenL)
				{
					p3=fourTimesLatticeMomentaComponentSquared[half_L]; //c2 - 32
					dummy_l2 += 32.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=0; // c3 -32
					dummy_l2 += 32.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
			if(evenL)
			{
				p2=fourTimesLatticeMomentaComponentSquared[half_L];
				{
					p3=fourTimesLatticeMomentaComponentSquared[half_L]; //d2 - 24
					dummy_l2 += 24.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=0; // b4 -48
					dummy_l2 += 48.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
			{
				p2=0;
				{
					p3=0; // d3 -24
					dummy_l2 += 24.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
		}
		if(evenL)
		{
			p1=fourTimesLatticeMomentaComponentSquared[half_L];
			{
				p2=fourTimesLatticeMomentaComponentSquared[half_L];
				{
					p3=fourTimesLatticeMomentaComponentSquared[half_L]; //c4 - 8
					dummy_l2 += 8.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=0; // b7 -24
					dummy_l2 += 24.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
			{
				p2=0;
				{
					p3=0; // b8 -24
					dummy_l2 += 24.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
		}
		{
			p1=0;
			{
				p2=0;
				{
					p3=0; // c5 - 8
					dummy_l2 += 8.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				dummy_l1+=dummy_l2;
			}
			dummy_l0+=dummy_l1;
		}
		sum+=dummy_l0;
	}
	dummy_l2=0.0;
	if(evenL)
	{
		p0=fourTimesLatticeMomentaComponentSquared[half_L];
		{
			p1=fourTimesLatticeMomentaComponentSquared[half_L];
			{
				p2=fourTimesLatticeMomentaComponentSquared[half_L];
				{
					p3=fourTimesLatticeMomentaComponentSquared[half_L]; //e2 - 1 
					dummy_l2 += 1.0/( p0 + p1 + p2 + p3 + massSquared);
				}
				{
					p3=0; // c7 -4
					dummy_l2 += 4.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
			{
				p2=0;
				{
					p3=0; // d4 -6
					dummy_l2 += 6.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
		}
		{
			p1=0;
			{
				p2=0;
				{
					p3=0; // c6 - 4
					dummy_l2 += 4.0/( p0 + p1 + p2 + p3 + massSquared);
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
					//ignore zero mode
// 					dummy_l2 += 1.0/( p0 + p1 + p2 + p3 + massSquared);
				}
			}
		}
		sum+=dummy_l2; //NOTE different dummy here
	}
	sum/=static_cast< double >(L0); sum/=static_cast< double >(L1); sum/=static_cast< double >(L2); sum/=static_cast< double >(L3);
	delete[] fourTimesLatticeMomentaComponentSquared;
	return sum;
}	
	
	
	
	
