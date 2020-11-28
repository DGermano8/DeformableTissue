#include "RandomRadialCellKiller3d.hpp"
#include "Debug.hpp"

/* Apoptosis for cells that are epithelial and lose contact with the basement membrane
 *
 */
RandomRadialCellKiller3d::RandomRadialCellKiller3d(AbstractCellPopulation<3>* pCrypt, double probabilityOfDeathInAnHour,
		double deathRadius, double centroidXCoordinate, double centroidYCoordinate)
    : AbstractCellKiller<3>(pCrypt),
      mProbabilityOfDeathInAnHour(probabilityOfDeathInAnHour),
      mDeathRadius(deathRadius),
      mCentroidXCoordinate(centroidXCoordinate),
      mCentroidYCoordinate(centroidYCoordinate)
{
}

RandomRadialCellKiller3d::~RandomRadialCellKiller3d()
{
}

void RandomRadialCellKiller3d::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
}

std::string RandomRadialCellKiller3d::GetOutputDirectory()
{
	return mOutputDirectory;
}

/** A method to return a vector that indicates which cells should be killed by anoikis
 * and which by compression-driven apoptosis
 */
std::vector<c_vector<unsigned,2> > RandomRadialCellKiller3d::RemoveByRandomSelection()
{
	MeshBasedCellPopulation<3>* p_tissue = static_cast<MeshBasedCellPopulation<3>*> (this->mpCellPopulation);
    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    std::vector<c_vector<unsigned,2> > cells_to_remove;
    c_vector<unsigned,2> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)
    bool outside_death_radius = false;				// Whether or not the cell is outside the death radius

	for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
	{
		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
		assert((!p_tissue->IsGhostNode(node_index)));

		// Initialise
		individual_node_information[0] = node_index;
		individual_node_information[1] = 0;

	    // We assume a constant time step
	    double death_prob_this_timestep = 1.0 - pow((1.0 - mProbabilityOfDeathInAnHour), SimulationTime::Instance()->GetTimeStep());
	    double x = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[0];
	    double y = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[1];

	    outside_death_radius = ( (x-mCentroidXCoordinate)*(x-mCentroidXCoordinate) + (y-mCentroidYCoordinate)*(y-mCentroidYCoordinate) > mDeathRadius );

	    if ((!cell_iter->HasApoptosisBegun() && (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false)
	    		&& (RandomNumberGenerator::Instance()->ranf() < death_prob_this_timestep) )
	    		&& (outside_death_radius) )
	    {
	    	cell_iter->StartApoptosis();
	    	individual_node_information[1] = 1;
	    }

		cells_to_remove.push_back(individual_node_information);
	}

	return cells_to_remove;
}

void RandomRadialCellKiller3d::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
    // We assume a constant time step
    double death_prob_this_timestep = 1.0 - pow((1.0 - mProbabilityOfDeathInAnHour), SimulationTime::Instance()->GetTimeStep());

    if (!pCell->HasApoptosisBegun() &&
        RandomNumberGenerator::Instance()->ranf() < death_prob_this_timestep)
    {
        pCell->StartApoptosis();
    }
}

/* Cell Killer that kills epithelial cells in the target zone (bounded by the specified x and y boundaries) by random selection
*/
void RandomRadialCellKiller3d::CheckAndLabelCellsForApoptosisOrDeath()
{
	MeshBasedCellPopulation<3>* p_tissue = static_cast<MeshBasedCellPopulation<3>*> (this->mpCellPopulation);
    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    std::vector<c_vector<unsigned,2> > cells_to_remove = RemoveByRandomSelection();

    // Loop over this vector and kill any cells that it tells you to
    for (unsigned i=0; i<cells_to_remove.size(); i++)
    {
    	if (cells_to_remove[i][1] == 1)
    	{
    		// Get cell associated to this node
    		CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
    		p_cell->Kill();
    	}
    }
}

void RandomRadialCellKiller3d::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ProbabilityOfDeathInAnHour>" << mProbabilityOfDeathInAnHour << "</ProbabilityOfDeathInAnHour> \n";
    *rParamsFile << "\t\t\t<DeathRadius>" << mDeathRadius << "</DeathRadius> \n";

    // Call direct parent class
    AbstractCellKiller<3>::OutputCellKillerParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RandomRadialCellKiller3d)
