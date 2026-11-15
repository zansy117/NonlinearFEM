#include "nonlinearfem/Assembly.h"

#include <stdexcept>
#include <cmath>
#include <limits>

namespace NonlinearFEM {
    namespace {
        void validateQ4MeshGeometry (
            const Q4Mesh& mesh,
            double thickness
        ) {
            const std::size_t N = mesh.nodes.size();

            if (N == 0 || mesh.elements.empty()) {
                throw std::invalid_argument(
                    "Q4 mesh must contain nodes and elements."
                );
            }

            if (N > std::numeric_limits<std::size_t>::max() / 2) {
                throw std::overflow_error(
                    "Global DOF count overflow."
                );
            }


            if (!std::isfinite(thickness) || thickness <= 0.0) {
                throw std::invalid_argument(
                    "Invalid reference thickness."
                );
            }

            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < 2; ++j) {
                    if (!std::isfinite(mesh.nodes[i][j])) {
                        throw std::invalid_argument(
                            "Non-finite nodal coordinate."
                        );
                    }


                }
            }

            for (const auto& element : mesh.elements) {
                for (std::size_t i = 0; i < kQ4NodeCount; ++i) {
                    if (element[i] >= N) {
                        throw std::out_of_range(
                            "Q4 element contains an invalid node index."
                        );
                    }

                    for (std::size_t j = i + 1;j < kQ4NodeCount; ++j) {
                        if (element[i] == element[j]) {
                            throw std::invalid_argument(
                                "Q4 element contains duplicate nodes."
                            );
                        }
                    }
                }
            }
        }

        void validateQ4Mesh (
            const Q4Mesh& mesh,
            const GlobalDisplacement& displacement,
            double thickness
        ) {
            const std::size_t N = mesh.nodes.size();

            if (N == 0 || mesh.elements.empty()) {
                throw std::invalid_argument(
                    "Q4 mesh must contain nodes and elements."
                );
            }

            if (N > std::numeric_limits<std::size_t>::max() / 2) {
                throw std::overflow_error(
                    "Global DOF count overflow."
                );
            }

            if (displacement.size() != 2 * N) {
                throw std::invalid_argument(
                    "Global displacement size does not match mesh DOFs."
                );
            }

            if (!std::isfinite(thickness) || thickness <= 0.0) {
                throw std::invalid_argument(
                    "Invalid reference thickness."
                );
            }

            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < 2; ++j) {
                    if (!std::isfinite(mesh.nodes[i][j])) {
                        throw std::invalid_argument(
                            "Non-finite nodal coordinate."
                        );
                    }

                    if (!std::isfinite(displacement[2 * i + j])) {
                        throw std::invalid_argument(
                            "Non-finite nodal displacement."
                        );
                    }
                }
            }

            for (const auto& element : mesh.elements) {
                for (std::size_t i = 0; i < kQ4NodeCount; ++i) {
                    if (element[i] >= N) {
                        throw std::out_of_range(
                            "Q4 element contains an invalid node index."
                        );
                    }

                    for (std::size_t j = i + 1;j < kQ4NodeCount; ++j) {
                        if (element[i] == element[j]) {
                            throw std::invalid_argument(
                                "Q4 element contains duplicate nodes."
                            );
                        }
                    }
                }
            }
        }    
    }





    [[nodiscard]]
    GlobalAssemblyResult assembleGlobalSystem (
        const Q4Mesh& mesh,
        const GlobalDisplacement& displacement,
        const MaterialParameters& material ,
        double thickness
    ) {
        validateQ4Mesh(mesh, displacement, thickness);

        const size_t nodeCount = mesh.nodes.size();
        const size_t dofCount = 2 * nodeCount;

        GlobalAssemblyResult result{};

        result.internalForce.resize(dofCount, 0.0);

        result.tangentStiffness.assign(dofCount, std::vector<double>(dofCount, 0.0));

        for(std::size_t i = 0;i < mesh.elements.size() ;++i) {
            Q4NodalValues elementCoordinates{};
            Q4NodalValues elementDisplacements{};

            Q4ElementConnectivity element = mesh.elements[i];

            std::array<std::size_t, kQ4DofCount> position{};

            for(std::size_t j = 0 ;j < kQ4NodeCount ;++j) {
                const size_t id = element[j];

                elementCoordinates[j] = mesh.nodes[id];

                elementDisplacements[j][0] = displacement[2 * id];
                elementDisplacements[j][1] = displacement[2 * id + 1];

                position[2 * j] = 2 * id;
                position[2 * j + 1] = 2 * id + 1;
            }

            const Q4InternalForce elementInternalForce = computeQ4InternalForce(elementCoordinates ,elementDisplacements ,material ,thickness);

            const Q4TangentStiffness elementTangentStiffness = computeQ4TangentStiffness(elementCoordinates ,elementDisplacements ,material,thickness);

            for(std::size_t j = 0 ;j < kQ4DofCount ;++j) {
                result.internalForce[position[j]] += elementInternalForce[j];
            }

            for(std::size_t j = 0;j < kQ4DofCount ;++j) {
                for(std::size_t k = 0;k < kQ4DofCount ;++k) {
                    result.tangentStiffness[position[j]][position[k]] += elementTangentStiffness[j][k];
                }
            } 
        }
        return result;
    }

    [[nodiscard]] 
    GlobalInternalForce assembleInternalForce (
        const Q4Mesh& mesh,
        const GlobalDisplacement& displacement,
        const MaterialParameters& material ,
        double thickness
    ) {
        validateQ4Mesh(mesh, displacement, thickness);

        const size_t N = mesh.nodes.size();

        GlobalInternalForce internalForce{};

        internalForce.resize(2 * N, 0.0);


        for(std::size_t i = 0;i < mesh.elements.size() ;++i) {
            Q4NodalValues elementCoordinates{};
            Q4NodalValues elementDisplacements{};

            Q4ElementConnectivity element = mesh.elements[i];

            std::array<std::size_t, kQ4DofCount> position{};

            for(std::size_t j = 0 ;j < kQ4NodeCount ;++j) {
                const size_t id = element[j];

                elementCoordinates[j] = mesh.nodes[id];

                elementDisplacements[j][0] = displacement[2 * id];
                elementDisplacements[j][1] = displacement[2 * id + 1];

                position[2 * j] = 2 * id;
                position[2 * j + 1] = 2 * id + 1;
            }

            const Q4InternalForce elementInternalForce = computeQ4InternalForce(elementCoordinates ,elementDisplacements ,material ,thickness);


            for(std::size_t j = 0 ;j < kQ4DofCount ;++j) {
                internalForce[position[j]] += elementInternalForce[j];
            }
        }
        return internalForce;
    }

    GlobalLumpedMass assembleGlobalLumpedMass (
        const Q4Mesh& mesh,
        double referenceDensity,
        double thickness
    ) {
        validateQ4MeshGeometry(mesh, thickness);

        GlobalLumpedMass mass{};

        const size_t N = mesh.nodes.size();

        mass.assign(2 * N,0.0);

        for(std::size_t i = 0;i < mesh.elements.size() ;++i) {
            Q4NodalValues elementCoordinates{};
            Q4ElementConnectivity element = mesh.elements[i];

            std::array<std::size_t, kQ4DofCount> position{};

            for(std::size_t j = 0;j < kQ4NodeCount;++j) {
                const size_t id = element[j];

                elementCoordinates[j] = mesh.nodes[id];

                position[2 * j] = 2 * id;
                position[2 * j + 1] = 2 * id + 1;
            }

            const Q4LumpedMass q4Mass = computeQ4LumpedMass(elementCoordinates,referenceDensity,thickness);

            for(std::size_t j = 0 ;j < kQ4DofCount ;++j) {
                mass[position[j]] += q4Mass[j];
            }
        }
        return mass;
    }

    double assembleGlobalEnergy (
        const Q4Mesh& mesh,
        const GlobalDisplacement& displacement,
        const MaterialParameters& material,
        double thickness 
    ) {
        validateQ4Mesh(mesh, displacement, thickness);

        const size_t N = mesh.elements.size();

        double energy = 0.0;

        for(std::size_t i = 0;i < N;++i) {
            Q4NodalValues elementCoordinates{};
            Q4NodalValues elementDisplacements{};

            Q4ElementConnectivity element = mesh.elements[i];


            for(std::size_t j = 0 ;j < kQ4NodeCount ;++j) {
                const size_t id = element[j];

                elementCoordinates[j] = mesh.nodes[id];

                elementDisplacements[j][0] = displacement[2 * id];
                elementDisplacements[j][1] = displacement[2 * id + 1];
            }

            energy += computeQ4StrainEnergy(elementCoordinates,elementDisplacements,material,thickness);
        }

        if (!std::isfinite(energy) || energy < 0.0) {
            throw std::invalid_argument(
                "Global strain energy must be finite and nonnegative."
            );
        }

        return energy;
    }
}
