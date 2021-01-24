/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef AbstractPositionAndPhaseBasedCellCycleModel_HPP_
#define AbstractPositionAndPhaseBasedCellCycleModel_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include <boost/serialization/base_object.hpp>

#include <vector>

#include "AbstractCellCycleModel.hpp"
#include "CellCyclePhases.hpp"
#include "SimulationTime.hpp"


/**
 * The AbstractPositionAndPhaseBasedCellCycleModel contains basic information to all phase based cell-cycle models.
 * It handles assignment of aspects of cell cycle phase.
 *
 */
class AbstractPositionAndPhaseBasedCellCycleModel : public AbstractCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mCurrentCellCyclePhase;
        archive & mG1Duration;
        archive & mMinimumGapDuration;
        archive & mStemCellG1Duration;
        archive & mTransitCellG1Duration;
        archive & mSDuration;
        archive & mG2Duration;
        archive & mMDuration;
    }

protected:

    /** The phase of the cell cycle that this model is in (specified in CellCyclePhases.hpp) */
    CellCyclePhase mCurrentCellCyclePhase;

    /**
     * How long the G1 phase lasts for.
     * Not necessarily a fixed value.
     */
    double mG1Duration;

    /**
     * Minimum possible duration of either of the gap phases (G1 or G2).
     * Has units of hours.
     *
     * Used to guarantee a strictly positive duration in cell-cycle models that
     * use normal random deviates for G1 or G2 phases.
     */
    double mMinimumGapDuration;

    /**
     * Duration of G1 phase for stem cells.
     * May be used as a mean duration for stochastic cell-cycle models.
     *
     */
    double mStemCellG1Duration;

    /**
     * Duration of G1 phase for transit cells.
     * May be used as a mean duration for stochastic cell-cycle models.
     */
    double mTransitCellG1Duration;

    /**
     * Duration of S phase for all cell types.
     */
    double mSDuration;

    /**
     * Duration of G2 phase for all cell types.
     */
    double mG2Duration;

    /**
     * Duration of M phase for all cell types.
     */
    double mMDuration;

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    AbstractPositionAndPhaseBasedCellCycleModel(const AbstractPositionAndPhaseBasedCellCycleModel& rModel);

public:

    /**
     * Default constructor - creates an AbstractPositionAndPhaseBasedCellCycleModel.
     */
    AbstractPositionAndPhaseBasedCellCycleModel();

    /**
     * Destructor.
     */
    virtual ~AbstractPositionAndPhaseBasedCellCycleModel();

    /**
     * See AbstractCellCycleModel::ResetForDivision()
     *
     * @return whether the cell is ready to divide (enter M phase).
     */
    virtual bool ReadyToDivide();

    /** See AbstractCellCycleModel::ResetForDivision() */
    virtual void ResetForDivision();

    /**
     * Set the phase the cell-cycle model is currently in. This method is called
     * from ReadyToDivide() just prior to deciding whether to divide the cell,
     * based on how far through the cell cycle it is, i.e. whether it has
     * completed M, G1, S and G2 phases.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void UpdateCellCyclePhase()=0;

    /**
     * @return the current cell cycle phase
     */
    CellCyclePhase GetCurrentCellCyclePhase() const;

    /**
     * @return the duration of the G1 phase of the cell cycle
     */
    virtual double GetG1Duration() const;

    /**
     * @return mStemCellG1Duration
     */
    double GetStemCellG1Duration() const;

    /**
     * @return mTransitCellG1Duration
     */
    double GetTransitCellG1Duration() const;

    /**
     * @return mSDuration + mG2Duration + mMDuration
     */
    double GetSG2MDuration() const;

    /**
     * @return the duration of the S phase of the cell cycle mSDuration
     */
    virtual double GetSDuration() const;

    /**
     * @return the duration of the G2 phase of the cell cycle mG2Duration
     */
    virtual double GetG2Duration() const;

    /**
     * @return the duration of the M phase of the cell cycle mMDuration
     */
    virtual double GetMDuration() const;

    /**
     * Set mStemCellG1Duration.
     *
     * @param stemCellG1Duration  the new value of mStemCellG1Duration
     */
    virtual void SetStemCellG1Duration(double stemCellG1Duration);

    /**
     * Set mTransitCellG1Duration.
     *
     * @param transitCellG1Duration  the new value of mTransitCellG1Duration
     */
    virtual void SetTransitCellG1Duration(double transitCellG1Duration);

    /**
     * Set mSDuration.
     *
     * @param sDuration  the new value of mSDuration
     */
    void SetSDuration(double sDuration);

    /**
     * Set mG2Duration.
     *
     * @param g2Duration  the new value of mG2Duration
     */
    void SetG2Duration(double g2Duration);

    /**
     * Set mMDuration.
     *
     * @param mDuration  the new value of mMDuration
     */
    void SetMDuration(double mDuration);

    /**
     * @return the typical cell cycle duration for a transit cell, in hours.
     * This method is overridden in some subclasses.
     */
    virtual double GetAverageTransitCellCycleTime();

    /**
     * @return the typical cell cycle duration for a stem cell, in hours.
     * This method is overridden in some subclasses.
     */
    virtual double GetAverageStemCellCycleTime();

    /**
     * @return mMinimumGapDuration
     */
    double GetMinimumGapDuration() const;

    /**
     * Set mMinimumGapDuration
     *
     * @param minimumGapDuration the new value of mMinimumGapDuration
     */
    void SetMinimumGapDuration(double minimumGapDuration);

    /**
     * Outputs cell cycle model parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile)=0;
};

CLASS_IS_ABSTRACT(AbstractPositionAndPhaseBasedCellCycleModel)

#endif /*AbstractPositionAndPhaseBasedCellCycleModel_HPP_*/
