#ifndef EMITTER
#define EMITTER

#include <iostream>
#include <vector>

#include <eigen3/Eigen/Core>

#include <ngl/Mat4.h>
#include <ngl/Camera.h>

#include "Particle.h"
#include "AlembicExport.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------
/// @file Emitter.h
/// @brief Class which contains and controls the particles making up a simulated object
/// @author Ina M. Sorensen
/// @version 1.0
/// @date 25.06.16
///
/// @todo
//------------------------------------------------------------------------------------------------------------------------------------------------------

class Emitter
{
  friend class Grid;

public:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Constructor
  //----------------------------------------------------------------------------------------------------------------------
  Emitter();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Destructor. Removes particle pointers
  //----------------------------------------------------------------------------------------------------------------------
  ~Emitter();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Generate particles. Can't be done in constructor as requires simulation constants to be read in first.
  //----------------------------------------------------------------------------------------------------------------------
  void createParticles(int _noParticles, const std::vector<Eigen::Vector3f> &_particlePositions, const std::vector<float> &_particleMass, const std::vector<float> &_particleTemperature, const std::vector<float> &_particlePhase);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set material constants
  //----------------------------------------------------------------------------------------------------------------------
  void setStrainConstants(float _lameMuConstant, float _lameLambdaConstant, float _compressionLim, float _stretchLim, float _hardnessCoefficient);
  void setTemperatureConstants(float _heatCapSolid, float _heatCapFluid, float _heatCondSolid, float _heatCondFluid, float _latentHeat, float _transitionTemp);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set collision object
  //----------------------------------------------------------------------------------------------------------------------
  void setCollisionObject(float _xMin, float _xMax, float _yMin, float _yMax, float _zMin, float _zMax);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set render parameters
  //----------------------------------------------------------------------------------------------------------------------
  void setRenderParameters(std::string _shaderName, float _particleRadius);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get compression limit
  //----------------------------------------------------------------------------------------------------------------------
  inline float getCompressionLimit() const {return m_compressionLimit;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get stretch limit
  //----------------------------------------------------------------------------------------------------------------------
  inline float getStretchLimit() const {return m_stretchLimit;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get hardness constant
  //----------------------------------------------------------------------------------------------------------------------
  inline float getHardnessCoefficient() const {return m_hardnessCoefficient;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get Lame Mu Constant
  //----------------------------------------------------------------------------------------------------------------------
  inline float getLameMuConstant() const {return m_lameMuConstant;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get Lame Lambda Constant
  //----------------------------------------------------------------------------------------------------------------------
  inline float getLameLambdaConstant() const {return m_lameLambdaConstant;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get number of particles
  //----------------------------------------------------------------------------------------------------------------------
  inline int getNoParticles() const {return m_noParticles;}
//  //----------------------------------------------------------------------------------------------------------------------
//  /// @brief Get list of particles
//  //----------------------------------------------------------------------------------------------------------------------
//  inline std::vector<Particle*>* getParticlesList() {return &m_particles;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get transition temperature
  //----------------------------------------------------------------------------------------------------------------------
  inline float getTransitionTemperature() const {return m_transitionTemperature;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get heat capacity for fluid
  //----------------------------------------------------------------------------------------------------------------------
  inline float getHeatCapacityFluid() const {return m_heatCapacityFluid;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get heat capacity for solid
  //----------------------------------------------------------------------------------------------------------------------
  inline float getHeatCapacitySolid() const {return m_heatCapacitySolid;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get latent heat
  //----------------------------------------------------------------------------------------------------------------------
  inline float getLatentHeat() const {return m_latentHeat;}

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Preset particles for first time step. Applies plasticity and makes corrections to all deformation gradient
  /// dependent variables accordingly
  //----------------------------------------------------------------------------------------------------------------------
  void presetParticles(float _velocityContribAlpha, float _tempContribBeta);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Update particles.
  //----------------------------------------------------------------------------------------------------------------------
  void updateParticles(float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Render particles
  /// @param [in] _modelMatrixCamera gives the scene transformations that also need to be applied to particle positions
  /// @param [in] _camera is used to get the view and projection matrices to give to the shader.
  //----------------------------------------------------------------------------------------------------------------------
  void renderParticles(ngl::Mat4 _modelMatrixCamera, ngl::Camera *_camera, float _ambientTemp, float _heatSourceTemp);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Export particles
  /// @param [in] _alembicExporter: std::unique pointer to alembic exporter
  //----------------------------------------------------------------------------------------------------------------------
//  void exportParticles(std::unique_ptr <AlembicExport> _alembicExporter);
  void exportParticles(AlembicExport* _alembicExporter);


protected:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Lame constant mu. Depends on the material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_lameMuConstant;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Lame constant lambda. Depends on the material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_lameLambdaConstant;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Hardness coefficient
  //----------------------------------------------------------------------------------------------------------------------
  float m_hardnessCoefficient;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Compression limit sets the compression value above which compression goes from elastic to plastic+elastic.
  /// Depends on the material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_compressionLimit;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Stretch limit sets the stretch value above which stretch goes from elastic to plastic+elastic.
  /// Depends on the material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_stretchLimit;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Heat capacity of solid. Depends on the material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_heatCapacitySolid;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Heat capacity of fluid. Depends on the material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_heatCapacityFluid;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Heat conductivity of solid. Depends on the material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_heatConductivitySolid;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Heat conductivity of fluid. Depends on the material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_heatConductivityFluid;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Latent heat is the heat it takes for solid-fluid and fluid-solid conversion to happen.
  /// Depends on the material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_latentHeat;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Transition temperature of material being simulated
  //----------------------------------------------------------------------------------------------------------------------
  float m_transitionTemperature;

private:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Number of particles contained by emitter
  //----------------------------------------------------------------------------------------------------------------------
  int m_noParticles;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief List of particles contained by emitter
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<Particle*> m_particles;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Collision object boundaries
  //----------------------------------------------------------------------------------------------------------------------
  float m_xMin;
  float m_xMax;
  float m_yMin;
  float m_yMax;
  float m_zMin;
  float m_zMax;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Shader name used to set which shader to use when rendering particles
  //----------------------------------------------------------------------------------------------------------------------
  std::string m_particleShaderName;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Used to set size of the rendered particles
  //----------------------------------------------------------------------------------------------------------------------
  float m_particleRadius;

};

#endif // EMITTER

