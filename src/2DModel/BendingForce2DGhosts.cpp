#include "BendingForce2DGhosts.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include <cmath>

#include "Debug.hpp"

template<unsigned DIM>
BendingForce2DGhosts<DIM>::BendingForce2DGhosts()
    : AbstractForce<DIM>(),
	  mBendingCoefficient(0.01),
      mExponentParameter(2.0),
      mDomainWidth(100.0),
      mSetNonZeroTargetCurvatureRegion(false),
      mNonZeroTargetCurvature(DOUBLE_UNSET),
      mRadiusOfNonZeroTargetCurvatureRegion(DOUBLE_UNSET),
      mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion(DOUBLE_UNSET)

{
}

template<unsigned DIM>
BendingForce2DGhosts<DIM>::~BendingForce2DGhosts()
{
}

template<unsigned DIM>
void BendingForce2DGhosts<DIM>::SetBendingCoefficient(double bendingCoefficient)
{
    mBendingCoefficient = bendingCoefficient;
}
template<unsigned DIM>
double BendingForce2DGhosts<DIM>::GetBendingCoefficient()
{
    return mBendingCoefficient;
}

template<unsigned DIM>
void BendingForce2DGhosts<DIM>::SetExponentParameter(double exponentParameter)
{
    mExponentParameter = exponentParameter;
}
template<unsigned DIM>
double BendingForce2DGhosts<DIM>::GetExponentParameter()
{
    return mExponentParameter;
}
template<unsigned DIM>
void BendingForce2DGhosts<DIM>::SetDomainWidth(double domainWidth)
{
    mDomainWidth= domainWidth;
}
template<unsigned DIM>
double BendingForce2DGhosts<DIM>::GetDomainWidth()
{
    return mDomainWidth;
}

template<unsigned DIM>
void BendingForce2DGhosts<DIM>::SetCircularNonZeroTargetCurvatureRegion(bool setNonZeroTargetCurvatureRegion, double nonZeroTargetCurvature,
		double radiusOfNonZeroTargetCurvatureRegion, double xCoordinateOfCentreOfNonZeroTargetCurvatureRegion)
{
	mSetNonZeroTargetCurvatureRegion = setNonZeroTargetCurvatureRegion;
	mNonZeroTargetCurvature = nonZeroTargetCurvature;
	mRadiusOfNonZeroTargetCurvatureRegion = radiusOfNonZeroTargetCurvatureRegion;
	mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion = xCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

template<unsigned DIM>
std::vector<c_vector<double, DIM> > BendingForce2DGhosts<DIM>::FitPlaneAndFindImage(AbstractCellPopulation<DIM>&  rCellPopulation, std::vector<unsigned> second_neighs, c_vector<double, DIM> cell_i_loc)
{
    std::vector<c_vector<double, DIM> > image_location_second_neighbours;



    double radius_of_curvature = 1.0/mNonZeroTargetCurvature;
    
    // This only works if we have normal other than (1,0), gotta do the inverse for other case...
    double x_bar = 0.0;
    double y_bar = 0.0;
    double xy_dot = 0.0;
    double xx_dot = 0.0;

    for (unsigned i=0; i<second_neighs.size(); i++)
    {
        unsigned cell_i_ext =second_neighs[i];
        CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(cell_i_ext);
        c_vector<double, DIM> epithelial_location = rCellPopulation.GetLocationOfCellCentre(p_cell_i_ext);
        
        x_bar += epithelial_location[0];
        y_bar += epithelial_location[1];
        
        xy_dot += epithelial_location[0]*epithelial_location[1];
        xx_dot += epithelial_location[0]*epithelial_location[0];

    }
    x_bar = 1.0/(second_neighs.size()) * x_bar;
    y_bar = 1.0/(second_neighs.size()) * y_bar;
    
    double fit_slope = (xy_dot - second_neighs.size()*x_bar*y_bar)/(xx_dot - second_neighs.size()*pow(x_bar,2));
    
    if(sqrt(pow(fit_slope,2)) < pow(10,-10))
    {
        fit_slope = 0.0;
    }
    double fit_int = cell_i_loc[1] - fit_slope*cell_i_loc[0];
    // PRINT_VARIABLE(fit_slope);
    if (fit_slope != 0.0)
    {
        double normal_slope = -1.0/fit_slope;
        double mult_slope;

        if (normal_slope < 0)
        {
            mult_slope = -1.0;
        }
        else
        {
            mult_slope = +1.0;
        }
        // PRINT_2_VARIABLES(mult_slope,normal_slope);

        double c0 = cell_i_loc[1] - normal_slope*cell_i_loc[0];

        double x0 = (cell_i_loc[0]*(1+pow(normal_slope,2)) + mult_slope*sqrt(pow(radius_of_curvature,2)*(1+pow(normal_slope,2))))/(1+pow(normal_slope,2));
        double y0 = normal_slope*x0 + c0;

        // PRINT_2_VARIABLES(x0,y0);
        
        for (unsigned i=0; i<second_neighs.size(); i++)
        {
            c_vector<double, DIM> temp_image = zero_vector<double>(DIM);

            unsigned cell_i_ext =second_neighs[i];
            CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(cell_i_ext);
            c_vector<double, DIM> epithelial_location = rCellPopulation.GetLocationOfCellCentre(p_cell_i_ext);

            double cij = epithelial_location[1] - normal_slope*epithelial_location[0];
            double xp = (cij-fit_int)/(fit_slope-normal_slope);
            double yp = normal_slope*xp + cij;

            double xc = (((y0-cij)*normal_slope+x0 - mult_slope*sqrt((1+pow(normal_slope,2))*pow(radius_of_curvature,2) - pow(cij+normal_slope*x0-y0,2)))/(1+pow(normal_slope,2)));
            double yc = normal_slope*xc +cij;

            temp_image[0] = epithelial_location[0] - (xc - xp);
            temp_image[1] = epithelial_location[1] - (yc - yp);

            // PRINT_VECTOR(temp_image)
            image_location_second_neighbours.push_back(temp_image);
        }

    }
    else
    {
        double x0 = cell_i_loc[0];
        double y0 = cell_i_loc[1] + radius_of_curvature;

        for (unsigned i=0; i<second_neighs.size(); i++)
        {
            c_vector<double, DIM> temp_image = zero_vector<double>(DIM);

            unsigned cell_i_ext =second_neighs[i];
            CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(cell_i_ext);
            c_vector<double, DIM> epithelial_location = rCellPopulation.GetLocationOfCellCentre(p_cell_i_ext);

            double xc = epithelial_location[0];
            double yc = y0 - sqrt(pow(radius_of_curvature,2) - pow(x0 - epithelial_location[0],2));

            double xp = epithelial_location[0];
            double yp = epithelial_location[1];

            temp_image[0] = epithelial_location[0];
            temp_image[1] = epithelial_location[1] - (yc - yp);

            // PRINT_VECTOR(temp_image)
            image_location_second_neighbours.push_back(temp_image);
        
        }
    }



    // c_vector<double, 4> ATA = zero_vector<double>(DIM*DIM);
    // c_vector<double, 2> ATz = zero_vector<double>(DIM);
    // c_vector<double, 4> iATA = zero_vector<double>(DIM*DIM);
    // c_vector<double, 2> normal_vector = zero_vector<double>(DIM);

    

    // for (unsigned i=0; i<second_neighs.size(); i++)
    // {
    //     unsigned cell_i_ext =second_neighs[i];
    //     CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(cell_i_ext);
    //     c_vector<double, DIM> epithelial_location = rCellPopulation.GetLocationOfCellCentre(p_cell_i_ext);
    //     // Fixed for periodicity
    //     // if(cell_i_loc[0] > 0.5*mDomainWidth)
    //     // {
    //     //     if(epithelial_location[0] < 0.25*mDomainWidth)
    //     //     {
    //     //         epithelial_location[0] = epithelial_location[0] + mDomainWidth;
    //     //     }
    //     // }
    //     // else
    //     // {
    //     //     if(epithelial_location[0] > 0.25*mDomainWidth)
    //     //     {
    //     //         epithelial_location[0] = epithelial_location[0] - mDomainWidth;
    //     //     }
    //     // }

    //     // Find Transpose(A)*A
    //     ATA[0] += epithelial_location[0]*epithelial_location[0];
    //     ATA[1] += epithelial_location[0];
    //     ATA[2] += epithelial_location[0];
    //     ATA[3] += 1;

    //     // Calculate Transpose(A)*z
    //     ATz[0] += epithelial_location[0]*epithelial_location[1];
    //     ATz[1] += epithelial_location[1];
    // }

    // // determinant of Transpose(A)*A
    // double ATA_det = ATA[0]*ATA[3] - ATA[2]*ATA[1];
    // if ( (ATA_det >= pow(10,-10)) || (ATA_det <= -1.0*pow(10,-10)) )
    // {
    //     // Calculate the inverse of Transpose(A)*A
    //     iATA[0] = (1.0/ATA_det)*ATA[3];
    //     iATA[1] = (1.0/ATA_det)*(-1.0*ATA[1]);
    //     iATA[2] = (1.0/ATA_det)*(-1.0*ATA[2]);
    //     iATA[3] = (1.0/ATA_det)*ATA[0];
                    
    //     // Calculate  normal = inverse(Transpose(A)*A)*(Transpose(A)*z)
    //     normal_vector[0] = iATA[0]*ATz[0] + iATA[1]*ATz[1];
    //     normal_vector[1] = iATA[2]*ATz[0] + iATA[3]*ATz[1];

    // }
    // normal_vector = normal_vector/(sqrt(pow(normal_vector[0],2)+pow(normal_vector[1],2)));
    
    // if( (0 - cell_i_loc[0])*normal_vector[0] + (10.0 - cell_i_loc[1])*normal_vector[1] < 0)
    // {
    //     normal_vector = -1.0*normal_vector;
    // }
    // PRINT_VECTOR(normal_vector);
    // // PRINT_VARIABLE(second_neighs.size());

    // // n = (0,1)
    // if(sqrt(pow(normal_vector[0],2))< 0.01 && sqrt(pow(normal_vector[1],2))> 1.0-0.01)
    // {
    //     double cirle_x = cell_i_loc[0] +radius_of_curvature*normal_vector[0];
    //     double cirle_y = cell_i_loc[1] +radius_of_curvature*normal_vector[1];
    //     // double slope_of_normal = normal_vector[1]/normal_vector[0];
                    
    //     for (unsigned i=0; i<second_neighs.size(); i++)
    //     {
    //         c_vector<double, DIM> temp_point_on_sphere = zero_vector<double>(DIM);
    //         c_vector<double, DIM> temp_image = zero_vector<double>(DIM);

    //         unsigned cell_i_ext =second_neighs[i];
    //         CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(cell_i_ext);
    //         c_vector<double, DIM> epithelial_location = rCellPopulation.GetLocationOfCellCentre(p_cell_i_ext);

    //         // if(cell_i_loc[0] > 0.5*mDomainWidth)
    //         // {
    //         //     if(epithelial_location[0] < 0.25*mDomainWidth)
    //         //     {
    //         //         epithelial_location[0] = epithelial_location[0] + mDomainWidth;
    //         //     }
    //         // }
    //         // else
    //         // {
    //         //     if(epithelial_location[0] > 0.25*mDomainWidth)
    //         //     {
    //         //         epithelial_location[0] = epithelial_location[0] - mDomainWidth;
    //         //     }
    //         // }
                      
    //         temp_point_on_sphere[0] = epithelial_location[0];
    //         temp_point_on_sphere[1] = cirle_y - sqrt(pow(radius_of_curvature,2) - pow(cirle_x - epithelial_location[0] ,2) );

    //         double distance_between_sphere_and_loc = sqrt(pow(epithelial_location[0]-temp_point_on_sphere[0],2) + pow(epithelial_location[1]-temp_point_on_sphere[1],2));

    //         temp_image[0] = epithelial_location[0] - distance_between_sphere_and_loc*normal_vector[0];
    //         temp_image[1] = epithelial_location[1] - distance_between_sphere_and_loc*normal_vector[1];

    //         image_location_second_neighbours.push_back(temp_image);
    //     }
    // }
    // // n = (1,0)
    // else if(sqrt(pow(normal_vector[1],2))< pow(10,-5) && sqrt(pow(normal_vector[0],2))> 1.0-pow(10,-5))
    // {
    //     // do the other thing!
    // }
    // else
    // {
    //     double cirle_x = cell_i_loc[0] +radius_of_curvature*normal_vector[0];
    //     double cirle_y = cell_i_loc[1] +radius_of_curvature*normal_vector[1];
    //     double slope_of_normal = normal_vector[1]/normal_vector[0];
                    
    //     for (unsigned i=0; i<second_neighs.size(); i++)
    //     {
    //         c_vector<double, DIM> temp_point_on_sphere = zero_vector<double>(DIM);
    //         c_vector<double, DIM> temp_image = zero_vector<double>(DIM);

    //         unsigned cell_i_ext =second_neighs[i];
    //         CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(cell_i_ext);
    //         c_vector<double, DIM> epithelial_location = rCellPopulation.GetLocationOfCellCentre(p_cell_i_ext);

    //         // if(cell_i_loc[0] > 0.5*mDomainWidth)
    //         // {
    //         //     if(epithelial_location[0] < 0.25*mDomainWidth)
    //         //     {
    //         //         epithelial_location[0] = epithelial_location[0] + mDomainWidth;
    //         //     }
    //         // }
    //         // else
    //         // {
    //         //     if(epithelial_location[0] > 0.25*mDomainWidth)
    //         //     {
    //         //         epithelial_location[0] = epithelial_location[0] - mDomainWidth;
    //         //     }
    //         // }

    //         double intersct_of_point = epithelial_location[1] - slope_of_normal*epithelial_location[0];
    //         double fit_a = pow(slope_of_normal,2)+1.0;
    //         double fit_b = 2*slope_of_normal*(intersct_of_point - cirle_y) - 2*cirle_x;
    //         double fit_c = pow(cirle_x,2) + pow(intersct_of_point - cirle_y,2) - pow(radius_of_curvature,2);
                        
    //         temp_point_on_sphere[0] = (-fit_b - sqrt(pow(fit_b,2) - 4*fit_a*fit_c) )/(2*fit_a);
    //         temp_point_on_sphere[1] = slope_of_normal*temp_point_on_sphere[0] + intersct_of_point;

    //         double distance_between_sphere_and_loc = sqrt(pow(epithelial_location[0]-temp_point_on_sphere[0],2) + pow(epithelial_location[1]-temp_point_on_sphere[1],2));

    //         temp_image[0] = epithelial_location[0] - distance_between_sphere_and_loc*normal_vector[0];
    //         temp_image[1] = epithelial_location[1] - distance_between_sphere_and_loc*normal_vector[1];

    //         image_location_second_neighbours.push_back(temp_image);
    //     }
    // }

    return image_location_second_neighbours;
}

template<unsigned DIM>
c_vector<double, DIM> BendingForce2DGhosts<DIM>::GetForceDueToDiscreteCurvature(AbstractCellPopulation<DIM>&  rCellPopulation, unsigned cell_i, std::vector<unsigned> first_neighs, std::vector<unsigned> second_neighs, std::vector<c_vector<double, DIM> > imag_loc)
{
    c_vector<double, DIM> force_due_to_curvature = zero_vector<double>(DIM);
    MeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

    // TRACE("curve hey 0");
    for(unsigned j=0; j<first_neighs.size(); j++)
	{
		double ang_sum = 0.0;
		c_vector<double, DIM> grad_i_phi_j_el = zero_vector<double>(DIM);;


		int cell_a = first_neighs[j];
		c_vector<double, DIM> a_loc, b_loc, c_loc;
        // TRACE("curve hey 1");
        // PRINT_VARIABLE(imag_loc.size());
        for(unsigned i=0; i<second_neighs.size(); i++)
        {
            if(second_neighs[i] == cell_a)
            {
                a_loc[0] = imag_loc[i][0];
                a_loc[1] = imag_loc[i][1];
            }
        } 
        // TRACE("curve hey 2");

        std::set<unsigned> all_first_neighbours = rCellPopulation.GetNeighbouringNodeIndices(cell_a);
        std::vector<unsigned> first_neighbours_epithelial;
        for (std::set<unsigned>::iterator iter = all_first_neighbours.begin();
        iter != all_first_neighbours.end();
        ++iter)
        {
            unsigned neighbour_node_global_index = *iter;
            if(!(p_tissue->IsGhostNode(neighbour_node_global_index)))
            {
                CellPtr p_cell_neigh = rCellPopulation.GetCellUsingLocationIndex(*iter);

                if(p_cell_neigh->GetMutationState()->IsType<WildTypeCellMutationState>())
                {
                    first_neighbours_epithelial.push_back(neighbour_node_global_index);
                }
            }
        }
        // PRINT_VARIABLE(first_neighbours_epithelial.size());
        int cell_b = first_neighbours_epithelial[0];
        int cell_c = first_neighbours_epithelial[1];
        

        for(unsigned i=0; i<second_neighs.size(); i++)
        {
            if(second_neighs[i] == cell_b)
            {
                b_loc[0] = imag_loc[i][0];
                b_loc[1] = imag_loc[i][1];
            }
            if(second_neighs[i] == cell_c)
            {
                c_loc[0] = imag_loc[i][0];
                c_loc[1] = imag_loc[i][1];
            }
        }
        


		// c_vector<double, DIM> vect_ab = rCellPopulation.rGetMesh().GetVectorFromAtoB(a_loc,b_loc);
		// c_vector<double, DIM> vect_ac = rCellPopulation.rGetMesh().GetVectorFromAtoB(a_loc,c_loc);
        c_vector<double, DIM> vect_ab = b_loc - a_loc;
		c_vector<double, DIM> vect_ac = c_loc - a_loc;

		// double mag_ab = sqrt(pow(b_loc[0] - a_loc[0],2) + pow(b_loc[1] - a_loc[1],2) );
		// double mag_ac = sqrt(pow(c_loc[0] - a_loc[0],2) + pow(c_loc[1] - a_loc[1],2) );
        double mag_ab = sqrt(pow(vect_ab[0],2) + pow(vect_ab[1],2) );
		double mag_ac = sqrt(pow(vect_ac[0],2) + pow(vect_ac[1],2) );

		double abDotac = (vect_ab[0])*(vect_ac[0]) + (vect_ab[1])*(vect_ac[1]);

		double cos_abc = (abDotac)/(mag_ab*mag_ac);

		if (!isnan(ang_sum))
		{
			if (cos_abc >= 0.99999999 || cos_abc == 1)
			{
				ang_sum = 0.0;
				cos_abc = 0.9999;
			}
			else if (cos_abc <= -0.99999999 || cos_abc == -1)
			{
				ang_sum = 3.141592653589793;
				cos_abc = -0.9999;
			}
			else
			{
                c_vector<double, DIM> d_loc;
                d_loc[0] = a_loc[0];
                d_loc[1] = a_loc[1];

                double mean_x = 1.0/3.0*(a_loc[0] + b_loc[0] + c_loc[0]);
                double var_x = 1.0/3.0*(pow(a_loc[0]-mean_x,2) + pow(b_loc[0]-mean_x,2) + pow(c_loc[0]-mean_x,2));

                double mean_y = 1.0/3.0*(a_loc[1] + b_loc[1] + c_loc[1]);
                double var_y = 1.0/3.0*(pow(a_loc[1]-mean_y,2) + pow(b_loc[1]-mean_y,2) + pow(c_loc[1]-mean_y,2));
                // std::cout<<"\n";
                if(var_y < 0.1)//small y variance
                {
                    // TRACE("1");
                    d_loc[0] = a_loc[0];
                    d_loc[1] = a_loc[1] - 20.0;

                }
                else if(var_x < 0.1)
                {
                    // TRACE("2");
                    d_loc[0] = a_loc[0] - 20.0;
                    d_loc[1] = a_loc[1];
                }
                else
                {
                    // TRACE("3");
                    d_loc[0] = a_loc[0] - 20.0;
                    d_loc[1] = a_loc[1] + 20.0;
                }
                double mag_ab_1 = sqrt(pow(b_loc[0] - a_loc[0],2) + pow(b_loc[1] - a_loc[1],2) );
                double mag_ad_1 = sqrt(pow(d_loc[0] - a_loc[0],2) + pow(d_loc[1] - a_loc[1],2) );
                double abDotad_1 = (b_loc[0] - a_loc[0])*(d_loc[0] - a_loc[0]) + (b_loc[1] - a_loc[1])*(d_loc[1] - a_loc[1]);
                double cos_abd_1 = (abDotad_1)/(mag_ab_1*mag_ad_1);

                double mag_ad_2 = sqrt(pow(d_loc[0] - a_loc[0],2) + pow(d_loc[1] - a_loc[1],2) );
	            double mag_ac_2 = sqrt(pow(c_loc[0] - a_loc[0],2) + pow(c_loc[1] - a_loc[1],2) );
		        double adDotac_2 = (d_loc[0] - a_loc[0])*(c_loc[0] - a_loc[0]) + (d_loc[1] - a_loc[1])*(c_loc[1] - a_loc[1]);
		        double cos_acd_2 = (adDotac_2)/(mag_ad_2*mag_ac_2);

			    ang_sum = acos(cos_abd_1) + acos(cos_acd_2);

			}
		}
		else if (isnan(ang_sum))
		{
			cos_abc = 0;
		}

		c_vector<double, DIM> grad_hold = zero_vector<double>(DIM);;

        // PRINT_VARIABLE(cos(ang_sum));

		// Calculate grad_i phi_i
        if(pow(cos(ang_sum),2) <= 0.999)
        {
            
            // PRINT_VECTOR(a_loc);
            // PRINT_VECTOR(b_loc);
            // PRINT_VECTOR(c_loc);
            // PRINT_VARIABLE(ang_sum);
            if(cell_a == cell_i)
            {
                // TRACE("a");
                grad_hold = ( -1.0/sqrt(1.0 - pow(cos(ang_sum),2)) )*(1.0/(mag_ab*mag_ac))*(-1.0*vect_ab - vect_ac + (abDotac)*( (vect_ab)/pow(mag_ab,2) +(vect_ac)/pow(mag_ac,2) ) );
                grad_i_phi_j_el += grad_hold;
            }
                    
            else if(cell_b == cell_i)
            {
                // TRACE("b");
                grad_hold = 1.0*(-1.0/(sqrt(1.0 - pow(cos(ang_sum),2))))*(1.0/(mag_ab*mag_ac))*( (vect_ac) - (vect_ab)*((abDotac))/pow(mag_ab,2) ) ;
                grad_i_phi_j_el += grad_hold;
            }
            else if(cell_c == cell_i)
            {
                // TRACE("c");
                grad_hold = 1.0*(-1.0/(sqrt(1.0 -  pow(cos(ang_sum),2))))*(1.0/(mag_ab*mag_ac))*( (vect_ab) - (vect_ac)*((abDotac))/pow(mag_ac,2) ) ;
                grad_i_phi_j_el += grad_hold;
            }

		double sign_of_anlge = (round((3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) >= 0) - (round((3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) < 0);
		c_vector<double, DIM> force_due_to_curvature_i = 1.0*(mExponentParameter)*sign_of_anlge*pow(sqrt(pow( (3.141592653589793 - ang_sum) ,2)),(mExponentParameter-1.0))*grad_i_phi_j_el;
        
		force_due_to_curvature +=  force_due_to_curvature_i;
        }
        // PRINT_VARIABLE(ang_sum);

	}


    return force_due_to_curvature;
}

template<unsigned DIM>
void BendingForce2DGhosts<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    MeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

    double dt = SimulationTime::Instance()->GetTimeStep();

    // Iterate over the nodes
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

		CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(real_node_index);

        

        if(p_cell->GetMutationState()->IsType<WildTypeCellMutationState>())
        {
            std::set<unsigned> all_first_neighbours = rCellPopulation.GetNeighbouringNodeIndices(real_node_index);

            std::vector<unsigned> first_neighbours_epithelial;
            for (std::set<unsigned>::iterator iter = all_first_neighbours.begin();
            iter != all_first_neighbours.end();
            ++iter)
            {
                unsigned neighbour_node_global_index = *iter;
                if(!(p_tissue->IsGhostNode(neighbour_node_global_index)))
                {
                    CellPtr p_cell_neigh = rCellPopulation.GetCellUsingLocationIndex(*iter);

                    if(p_cell_neigh->GetMutationState()->IsType<WildTypeCellMutationState>())
                    {
                        first_neighbours_epithelial.push_back(neighbour_node_global_index);
                    }
                }
            }
            first_neighbours_epithelial.push_back(real_node_index);
            std::vector<unsigned> second_neighbours_epithelial;

            for(unsigned j=0; j<first_neighbours_epithelial.size(); j++)
			{
				unsigned neigh_j_index = first_neighbours_epithelial[j];

				std::set<unsigned> all_first_neighbours_j = rCellPopulation.GetNeighbouringNodeIndices(neigh_j_index);
                for (std::set<unsigned>::iterator iter = all_first_neighbours_j.begin();
                iter != all_first_neighbours_j.end();
                ++iter)
                {
                    unsigned neighbour_node_global_index = *iter;
                    if(!(p_tissue->IsGhostNode(neighbour_node_global_index)))
                    {
                        CellPtr p_cell_neigh = rCellPopulation.GetCellUsingLocationIndex(*iter);

                        if(p_cell_neigh->GetMutationState()->IsType<WildTypeCellMutationState>())
                        {
                            second_neighbours_epithelial.push_back(neighbour_node_global_index);
                        }
                    }
                }

			}
            std::sort(second_neighbours_epithelial.begin(), second_neighbours_epithelial.end());
			second_neighbours_epithelial.erase(std::unique(second_neighbours_epithelial.begin(), second_neighbours_epithelial.end()), second_neighbours_epithelial.end());

            // PRINT_VARIABLE(second_neighbours_epithelial.size());

            c_vector<double, DIM> real_node_location_i = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            // TRACE("HEY 0");
            // Here we find the image of the nodes:
            std::vector<c_vector<double, DIM> > image_location_second_neighbours;
			bool epithelial_cell_in_non_target_curvature_region = ( sqrt(pow((real_node_location_i[0] - mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion),2)) < mRadiusOfNonZeroTargetCurvatureRegion);

            double beta_paramater = mBendingCoefficient;
            if(epithelial_cell_in_non_target_curvature_region)
            {
                // PRINT_VARIABLE(real_node_location_i[0]);
                // beta_paramater = -1.0*beta_paramater;
                image_location_second_neighbours = FitPlaneAndFindImage(rCellPopulation, second_neighbours_epithelial, real_node_location_i);
            }
            else
            {
                for(int neigh_j=0; neigh_j<second_neighbours_epithelial.size(); neigh_j++)
                {
                    int neigh_j_index = second_neighbours_epithelial[neigh_j];
                    CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(neigh_j_index);
                    c_vector<double, DIM> node_location_j = rCellPopulation.GetLocationOfCellCentre(p_cell_i_ext);
                    image_location_second_neighbours.push_back(node_location_j);
                }
            }
            c_vector<double, DIM> force_contribution= zero_vector<double>(DIM);
            if(image_location_second_neighbours.size() > 1)
            {
                if(real_node_location_i[0] < 2.5 || real_node_location_i[0] > mDomainWidth - 2.5)
                {
                    // do nothing
                    force_contribution[0] = 0.0;
                    force_contribution[1] = 0.0;
                }
                else
                {
                    force_contribution = GetForceDueToDiscreteCurvature(rCellPopulation, real_node_index, first_neighbours_epithelial, second_neighbours_epithelial, image_location_second_neighbours);
                }

            }
            if(isnan(force_contribution[0]) || isnan(force_contribution[1]))
            {
                force_contribution[0] = 0.0;
                force_contribution[1] = 0.0;
            }

            force_contribution = beta_paramater*force_contribution;

            rCellPopulation.GetNode(real_node_index)->AddAppliedForceContribution(force_contribution);
            // if(epithelial_cell_in_non_target_curvature_region)
            // {

            //     PRINT_VARIABLE(p_cell->GetCellId());
            //     PRINT_VECTOR(force_contribution);
            // }
            
        }
    }
    // std::cout<<"_____________________________________________\n\n";
}

template<unsigned DIM>
void BendingForce2DGhosts<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<BendingCoefficient>" << mBendingCoefficient << "</BendingCoefficient> \n";
    *rParamsFile << "\t\t\t<ExponentParameter>" << mExponentParameter << "</ExponentParameter> \n";
    *rParamsFile << "\t\t\t<DomainWidth>" << mDomainWidth << "</DomainWidth> \n";
    *rParamsFile << "\t\t\t<SetNonZeroTargetCurvatureRegion>" << mSetNonZeroTargetCurvatureRegion << "</SetNonZeroTargetCurvatureRegion> \n";
    *rParamsFile << "\t\t\t<NonZeroTargetCurvature>" << mNonZeroTargetCurvature << "</NonZeroTargetCurvature> \n";
    *rParamsFile << "\t\t\t<RadiusOfNonZeroTargetCurvatureRegion>" << mRadiusOfNonZeroTargetCurvatureRegion << "</RadiusOfNonZeroTargetCurvatureRegion> \n";
    *rParamsFile << "\t\t\t<XCoordinateOfCentreOfNonZeroTargetCurvatureRegion>" << mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion << "</XCoordinateOfCentreOfNonZeroTargetCurvatureRegion> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class BendingForce2DGhosts<1>;
template class BendingForce2DGhosts<2>;
template class BendingForce2DGhosts<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BendingForce2DGhosts)
