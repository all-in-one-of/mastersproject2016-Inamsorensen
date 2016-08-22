#include "Grid.h"

void Grid::updateParticleFromGrid_New(Emitter* _emitter, float _velocityContribAlpha, float _tempContribBeta)
{
  /* Outline
  ---------------------------------------------------------------------------------------------------------------------
  Loop over all particles

    Find the cells in the grid that the particle will contribute to
      For simplicity loop over the same number of cells as for cubic B Spline - This should be optimised

    Loop over cells
      Get quadratic stencil and its differential
      Get velocity and previous velocity
      Get temperature and previous temp

      Calc PIC velocity for each face
      Calc FLIP velocity for each face
      Add to particle

      Calc velocity gradient contribution for each face
      Add to particle

      Calc PIC temperature for cell centre
      Calc FLIP temperature for cell centre
      Add to particle

  ---------------------------------------------------------------------------------------------------------------------
  */

  //Set e_{a(i)} vectors
  Eigen::Vector3f e_x(1.0, 0.0, 0.0);
  Eigen::Vector3f e_y(0.0, 1.0, 0.0);
  Eigen::Vector3f e_z(0.0, 0.0, 1.0);

  //Get total number of particles
  int noParticles=_emitter->m_noParticles;

  //To calc position of particle, need origin of grid corner, not centre of first grid cell.
  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f gridEdgePosition=m_origin;
  gridEdgePosition(0)-=halfCellSize;
  gridEdgePosition(1)-=halfCellSize;
  gridEdgePosition(2)-=halfCellSize;

//  #pragma omp parallel for
  for (int particleItr=0; particleItr<noParticles; particleItr++)
  {
    Particle* particlePtr=_emitter->m_particles[particleItr];
    Eigen::Vector3f particlePosition=particlePtr->getPosition();
    Eigen::Vector3i particleIndex=MathFunctions::getParticleGridCell(particlePosition, m_cellSize, gridEdgePosition);

    //Loop over i+-2, j+-2, k+-2 to get cells that particle will contribute to
    int iParticle=particleIndex(0);
    int jParticle=particleIndex(1);
    int kParticle=particleIndex(2);

    for (int kIndex=(kParticle-2); kIndex<(kParticle+4); kIndex++)
    {
      for (int jIndex=(jParticle-2); jIndex<(jParticle+4); jIndex++)
      {
        for (int iIndex=(iParticle-2); iIndex<(iParticle+4); iIndex++)
        {
          //Check that not in outer cells or outside grid
          if (iIndex>=0 && iIndex<=(m_noCells-1) && jIndex>=0 && jIndex<=(m_noCells-1) && kIndex>=0 && kIndex<=(m_noCells-1))
          {
            //Get cell index
            int cellIndex=MathFunctions::getVectorIndex(iIndex, jIndex, kIndex, m_noCells);

            //Get tight quadratic stencil weight
            float weightCentre=0.0;
            float weightX=0.0;
            float weightY=0.0;
            float weightZ=0.0;
            calcWeight_tightQuadraticStencil(particlePosition, iIndex, jIndex, kIndex,
                                             weightCentre, weightX, weightY, weightZ);

            //Get tight quadratic stencil weight differentiated
            Eigen::Vector3f weightX_Diff;
            Eigen::Vector3f weightY_Diff;
            Eigen::Vector3f weightZ_Diff;
            weightX_Diff.setZero();
            weightY_Diff.setZero();
            weightZ_Diff.setZero();
            calcWeight_tightQuadraticStencil_Diff(particlePosition, iIndex, jIndex, kIndex,
                                                  weightX_Diff, weightY_Diff, weightZ_Diff);

            //Get velocity and previous velocity of faces
            float velocity_FaceX=m_cellFacesX[cellIndex]->m_velocity;
            float velocity_FaceY=m_cellFacesY[cellIndex]->m_velocity;
            float velocity_FaceZ=m_cellFacesZ[cellIndex]->m_velocity;
            float prevVelocity_FaceX=m_cellFacesX[cellIndex]->m_previousVelocity;
            float prevVelocity_FaceY=m_cellFacesY[cellIndex]->m_previousVelocity;
            float prevVelocity_FaceZ=m_cellFacesZ[cellIndex]->m_previousVelocity;

            //Get temperature and previous temperature
            float temperature=m_cellCentres[cellIndex]->m_temperature;
            float prevTemperature=m_cellCentres[cellIndex]->m_previousTemperature;


            //Only update if stencil is non-zero
            if (weightX!=0)
            {
              //Update velocity
              //PIC velocity
              float velocityPIC_FaceX=velocity_FaceX*weightX;

              //FLIP velocity
              float velocityFLIP_FaceX=(velocity_FaceX-prevVelocity_FaceX)*weightX;

              //Velocity contribution
              float velocityContribution_FaceX=(_velocityContribAlpha*velocityFLIP_FaceX)+((1.0-_velocityContribAlpha)*velocityPIC_FaceX);

              Eigen::Vector3f velContribVector=velocityContribution_FaceX*e_x;

              //Set up velocity gradient contribution
              Eigen::Matrix3f velGradContribution;
              velGradContribution.setZero();
              velGradContribution(0,0)=velocity_FaceX*weightX_Diff(0);
              velGradContribution(0,1)=velocity_FaceX*weightX_Diff(1);
              velGradContribution(0,2)=velocity_FaceX*weightX_Diff(2);

              //Update particle
              particlePtr->addParticleVelocity(velContribVector);
              particlePtr->addParticleVelocityGradient(velGradContribution);

            }

            if (weightY!=0)
            {
              //Update velocity
              //PIC velocity
              float velocityPIC_FaceY=velocity_FaceY*weightY;

              //FLIP velocity
              float velocityFLIP_FaceY=(velocity_FaceY-prevVelocity_FaceY)*weightY;

              //Velocity contribution
              float velocityContribution_FaceY=(_velocityContribAlpha*velocityFLIP_FaceY)+((1.0-_velocityContribAlpha)*velocityPIC_FaceY);

              Eigen::Vector3f velContribVector=velocityContribution_FaceY*e_y;

              //Set up velocity gradient contribution
              Eigen::Matrix3f velGradContribution;
              velGradContribution.setZero();
              velGradContribution(1,0)=velocity_FaceY*weightY_Diff(0);
              velGradContribution(1,1)=velocity_FaceY*weightY_Diff(1);
              velGradContribution(1,2)=velocity_FaceY*weightY_Diff(2);

              //Update particle
              particlePtr->addParticleVelocity(velContribVector);
              particlePtr->addParticleVelocityGradient(velGradContribution);
            }

            if (weightZ!=0)
            {
              //Update velocity
              //PIC velocity
              float velocityPIC_FaceZ=velocity_FaceZ*weightZ;

              //FLIP velocity
              float velocityFLIP_FaceZ=(velocity_FaceZ-prevVelocity_FaceZ)*weightZ;

              //Velocity contribution
              float velocityContribution_FaceZ=(_velocityContribAlpha*velocityFLIP_FaceZ)+((1.0-_velocityContribAlpha)*velocityPIC_FaceZ);

              Eigen::Vector3f velContribVector=velocityContribution_FaceZ*e_z;

              //Set up velocity gradient contribution
              Eigen::Matrix3f velGradContribution;
              velGradContribution.setZero();
              velGradContribution(2,0)=velocity_FaceZ*weightZ_Diff(0);
              velGradContribution(2,1)=velocity_FaceZ*weightZ_Diff(1);
              velGradContribution(2,2)=velocity_FaceZ*weightZ_Diff(2);

              //Update particle
              particlePtr->addParticleVelocity(velContribVector);
              particlePtr->addParticleVelocityGradient(velGradContribution);
            }

            if (weightCentre!=0)
            {
              //Temperature contribution
              //PIC temperature
              float temperaturePIC=temperature*weightCentre;

              //FLIP temperature
              float temperatureFLIP=(temperature-prevTemperature)*weightCentre;

              //Calculate temperature contribution
              float temperatureContribution=(_tempContribBeta*temperatureFLIP)+((1.0-_tempContribBeta)*temperaturePIC);

              //Update particle
              particlePtr->addParticleTemperature(temperatureContribution);
            }

          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcWeight_tightQuadraticStencil(Eigen::Vector3f _particlePosition, int _iIndex, int _jIndex, int _kIndex, float &o_weightCentre, float &o_weightFaceX, float &o_weightFaceY, float &o_weightFaceZ)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calculate interpolation weights for the particle and the cell with given _i, _j,_k

  Calculate i*, ix, iy, iz position vectors for the cell in question

  Calculate x=xp-xia values for these

  Pass in x and get interpolation weights:
       tightQuadraticStencil

  Return results
  ------------------------------------------------------------------------------------------------------
  */

  //Initialise weights to zero
  o_weightFaceX=0.0;
  o_weightFaceY=0.0;
  o_weightFaceZ=0.0;

  //Cell position
  float xPos=(_iIndex*m_cellSize)+m_origin(0);
  float yPos=(_jIndex*m_cellSize)+m_origin(1);
  float zPos=(_kIndex*m_cellSize)+m_origin(2);

  //Position vectors for centre and faces
  Eigen::Vector3f centreVector(xPos, yPos, zPos);

  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f faceXVector(xPos-halfCellSize, yPos, zPos);
  Eigen::Vector3f faceYVector(xPos, yPos-halfCellSize, zPos);
  Eigen::Vector3f faceZVector(xPos, yPos, zPos-halfCellSize);


  //Calculate posDifference for each face and cell centre
  Eigen::Vector3f centrePosDiff=_particlePosition-centreVector;
  Eigen::Vector3f faceXPosDiff=_particlePosition-faceXVector;
  Eigen::Vector3f faceYPosDiff=_particlePosition-faceYVector;
  Eigen::Vector3f faceZPosDiff=_particlePosition-faceZVector;

  //Centre
  float NxCentre_tightQS=MathFunctions::calcTightQuadraticStencil(centrePosDiff(0)/m_cellSize);
  float NyCentre_tightQS=MathFunctions::calcTightQuadraticStencil(centrePosDiff(1)/m_cellSize);
  float NzCentre_tightQS=MathFunctions::calcTightQuadraticStencil(centrePosDiff(2)/m_cellSize);

  //FaceX
  float NxFaceX_tightQS=MathFunctions::calcTightQuadraticStencil(faceXPosDiff(0)/m_cellSize);
  float NyFaceX_tightQS=MathFunctions::calcTightQuadraticStencil(faceXPosDiff(1)/m_cellSize);
  float NzFaceX_tightQS=MathFunctions::calcTightQuadraticStencil(faceXPosDiff(2)/m_cellSize);

  //FaceY
  float NxFaceY_tightQS=MathFunctions::calcTightQuadraticStencil(faceYPosDiff(0)/m_cellSize);
  float NyFaceY_tightQS=MathFunctions::calcTightQuadraticStencil(faceYPosDiff(1)/m_cellSize);
  float NzFaceY_tightQS=MathFunctions::calcTightQuadraticStencil(faceYPosDiff(2)/m_cellSize);

  //FaceZ
  float NxFaceZ_tightQS=MathFunctions::calcTightQuadraticStencil(faceZPosDiff(0)/m_cellSize);
  float NyFaceZ_tightQS=MathFunctions::calcTightQuadraticStencil(faceZPosDiff(1)/m_cellSize);
  float NzFaceZ_tightQS=MathFunctions::calcTightQuadraticStencil(faceZPosDiff(2)/m_cellSize);

  //Return interpolation weights
  o_weightCentre=NxCentre_tightQS*NyCentre_tightQS*NzCentre_tightQS;
  o_weightFaceX=NxFaceX_tightQS*NyFaceX_tightQS*NzFaceX_tightQS;
  o_weightFaceY=NxFaceY_tightQS*NyFaceY_tightQS*NzFaceY_tightQS;
  o_weightFaceZ=NxFaceZ_tightQS*NyFaceZ_tightQS*NzFaceZ_tightQS;

}

//----------------------------------------------------------------------------------------------------------------------

void Grid::calcWeight_tightQuadraticStencil_Diff(Eigen::Vector3f _particlePosition, int _iIndex, int _jIndex, int _kIndex, Eigen::Vector3f &o_weightFaceX, Eigen::Vector3f &o_weightFaceY, Eigen::Vector3f &o_weightFaceZ)
{
  /* Outline
  ------------------------------------------------------------------------------------------------------
  Calculate interpolation weights for the particle and the cell with given _i, _j,_k

  Calculate i*, ix, iy, iz position vectors for the cell in question

  Calculate x=xp-xia values for these

  Pass in x and get interpolation weights:
       tightQuadraticStencil_Diff

  Return results
  ------------------------------------------------------------------------------------------------------
  */

  //Initialise weights to zero
  o_weightFaceX.setZero();
  o_weightFaceY.setZero();
  o_weightFaceZ.setZero();

  //Cell position
  float xPos=(_iIndex*m_cellSize)+m_origin(0);
  float yPos=(_jIndex*m_cellSize)+m_origin(1);
  float zPos=(_kIndex*m_cellSize)+m_origin(2);

  //Position vectors for faces
  float halfCellSize=m_cellSize/2.0;
  Eigen::Vector3f faceXVector(xPos-halfCellSize, yPos, zPos);
  Eigen::Vector3f faceYVector(xPos, yPos-halfCellSize, zPos);
  Eigen::Vector3f faceZVector(xPos, yPos, zPos-halfCellSize);

  //Calculate posDifference for each face
  Eigen::Vector3f faceXPosDiff=_particlePosition-faceXVector;
  Eigen::Vector3f faceYPosDiff=_particlePosition-faceYVector;
  Eigen::Vector3f faceZPosDiff=_particlePosition-faceZVector;

  //Need to get the non-differentiated values as well
  float Nx_tightQS_FaceX=MathFunctions::calcTightQuadraticStencil(faceXPosDiff(0)/m_cellSize);
  float Ny_tightQS_FaceX=MathFunctions::calcTightQuadraticStencil(faceXPosDiff(1)/m_cellSize);
  float Nz_tightQS_FaceX=MathFunctions::calcTightQuadraticStencil(faceXPosDiff(2)/m_cellSize);

  //Face X
  float dNx_tightQS_FaceX=MathFunctions::calcTightQuadraticStencil_Diff(faceXPosDiff(0)/m_cellSize);
  float dNy_tightQS_FaceX=MathFunctions::calcTightQuadraticStencil_Diff(faceXPosDiff(1)/m_cellSize);
  float dNz_tightQS_FaceX=MathFunctions::calcTightQuadraticStencil_Diff(faceXPosDiff(2)/m_cellSize);

  o_weightFaceX(0)=dNx_tightQS_FaceX*Ny_tightQS_FaceX*Nz_tightQS_FaceX;
  o_weightFaceX(1)=dNy_tightQS_FaceX*Nx_tightQS_FaceX*Nz_tightQS_FaceX;
  o_weightFaceX(2)=dNz_tightQS_FaceX*Nx_tightQS_FaceX*Ny_tightQS_FaceX;

  o_weightFaceX*=(1.0/m_cellSize);


  //Face Y
  float Nx_tightQS_FaceY=MathFunctions::calcTightQuadraticStencil(faceYPosDiff(0)/m_cellSize);
  float Ny_tightQS_FaceY=MathFunctions::calcTightQuadraticStencil(faceYPosDiff(1)/m_cellSize);
  float Nz_tightQS_FaceY=MathFunctions::calcTightQuadraticStencil(faceYPosDiff(2)/m_cellSize);

  float dNx_tightQS_FaceY=MathFunctions::calcTightQuadraticStencil_Diff(faceYPosDiff(0)/m_cellSize);
  float dNy_tightQS_FaceY=MathFunctions::calcTightQuadraticStencil_Diff(faceYPosDiff(1)/m_cellSize);
  float dNz_tightQS_FaceY=MathFunctions::calcTightQuadraticStencil_Diff(faceYPosDiff(2)/m_cellSize);

  o_weightFaceY(0)=dNx_tightQS_FaceY*Ny_tightQS_FaceY*Nz_tightQS_FaceY;
  o_weightFaceY(1)=dNy_tightQS_FaceY*Nx_tightQS_FaceY*Nz_tightQS_FaceY;
  o_weightFaceY(2)=dNz_tightQS_FaceY*Nx_tightQS_FaceY*Ny_tightQS_FaceY;

  o_weightFaceY*=(1.0/m_cellSize);

  //Face Z
  float Nx_tightQS_FaceZ=MathFunctions::calcTightQuadraticStencil(faceZPosDiff(0)/m_cellSize);
  float Ny_tightQS_FaceZ=MathFunctions::calcTightQuadraticStencil(faceZPosDiff(1)/m_cellSize);
  float Nz_tightQS_FaceZ=MathFunctions::calcTightQuadraticStencil(faceZPosDiff(2)/m_cellSize);

  float dNx_tightQS_FaceZ=MathFunctions::calcTightQuadraticStencil_Diff(faceZPosDiff(0)/m_cellSize);
  float dNy_tightQS_FaceZ=MathFunctions::calcTightQuadraticStencil_Diff(faceZPosDiff(1)/m_cellSize);
  float dNz_tightQS_FaceZ=MathFunctions::calcTightQuadraticStencil_Diff(faceZPosDiff(2)/m_cellSize);

  o_weightFaceZ(0)=dNx_tightQS_FaceZ*Ny_tightQS_FaceZ*Nz_tightQS_FaceZ;
  o_weightFaceZ(1)=dNy_tightQS_FaceZ*Nx_tightQS_FaceZ*Nz_tightQS_FaceZ;
  o_weightFaceZ(2)=dNz_tightQS_FaceZ*Nx_tightQS_FaceZ*Ny_tightQS_FaceZ;

  o_weightFaceZ*=(1.0/m_cellSize);
}

//----------------------------------------------------------------------------------------------------------------------
