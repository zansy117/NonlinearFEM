#include "nonlinearfem/Assembly_3D.h"

#include <stdexcept>
#include <cmath>
#include <limits>


namespace NonlinearFEM {
    namespace {
        void validateH8MeshGeometry (
            const H8Mesh& mesh
        ) {
            const std::size_t N = mesh.nodes.size();

            if (N == 0 || mesh.elements.empty()) {
                throw std::invalid_argument(
                    "H8 mesh must contain nodes and elements."
                );
            }

            if (N > std::numeric_limits<std::size_t>::max() / 3) {
                throw std::overflow_error(
                    "Global DOF count overflow."
                );
            }

            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    if (!std::isfinite(mesh.nodes[i][j])) {
                        throw std::invalid_argument(
                            "Non-finite nodal coordinate."
                        );
                    }


                }
            }

            for (const auto& element : mesh.elements) {
                for (std::size_t i = 0; i < kH8NodeCount; ++i) {
                    if (element[i] >= N) {
                        throw std::out_of_range(
                            "H8 element contains an invalid node index."
                        );
                    }

                    for (std::size_t j = i + 1;j < kH8NodeCount; ++j) {
                        if (element[i] == element[j]) {
                            throw std::invalid_argument(
                                "H8 element contains duplicate nodes."
                            );
                        }
                    }
                }
            }
        }

        void validateH8Mesh (
            const H8Mesh& mesh,
            const GlobalDisplacement_3D& displacement
        ) {
            const std::size_t N = mesh.nodes.size();

            if (N == 0 || mesh.elements.empty()) {
                throw std::invalid_argument(
                    "H8 mesh must contain nodes and elements."
                );
            }

            if (N > std::numeric_limits<std::size_t>::max() / 3) {
                throw std::overflow_error(
                    "Global DOF count overflow."
                );
            }

            if (displacement.size() != 3 * N) {
                throw std::invalid_argument(
                    "Global displacement size does not match mesh DOFs."
                );
            }

            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    if (!std::isfinite(mesh.nodes[i][j])) {
                        throw std::invalid_argument(
                            "Non-finite nodal coordinate."
                        );
                    }

                    if (!std::isfinite(displacement[3 * i + j])) {
                        throw std::invalid_argument(
                            "Non-finite nodal displacement."
                        );
                    }
                }
            }

            for (const auto& element : mesh.elements) {
                for (std::size_t i = 0; i < kH8NodeCount; ++i) {
                    if (element[i] >= N) {
                        throw std::out_of_range(
                            "H8 element contains an invalid node index."
                        );
                    }

                    for (std::size_t j = i + 1;j < kH8NodeCount; ++j) {
                        if (element[i] == element[j]) {
                            throw std::invalid_argument(
                                "H8 element contains duplicate nodes."
                            );
                        }
                    }
                }
            }
        }    
    }

    

    GlobalAssemblyResult_3D assembleGlobalSystem_3D (
        const H8Mesh& mesh,
        const GlobalDisplacement_3D& displacement,
        const MaterialParameters& material
    ) {
        validateH8Mesh(mesh, displacement);

        const size_t nodeCount = mesh.nodes.size();
        const size_t dofCount = 3 * nodeCount;

        GlobalAssemblyResult_3D result{};

        result.internalForce.resize(dofCount, 0.0);

        result.tangentStiffness.assign(dofCount, std::vector<double>(dofCount, 0.0));

        for(std::size_t i = 0;i < mesh.elements.size() ;++i) {
            H8NodalValues elementCoordinates{};
            H8NodalValues elementDisplacements{};

            H8ElementConnectivity element = mesh.elements[i];

            std::array<std::size_t, kH8DofCount> position{};

            for(std::size_t j = 0 ;j < kH8NodeCount ;++j) {
                const size_t id = element[j];

                elementCoordinates[j] = mesh.nodes[id];

                elementDisplacements[j][0] = displacement[3 * id];
                elementDisplacements[j][1] = displacement[3 * id + 1];
                elementDisplacements[j][2] = displacement[3 * id + 2];

                position[3 * j] = 3 * id;
                position[3 * j + 1] = 3 * id + 1;
                position[3 * j + 2] = 3 * id + 2;
            }

            const H8InternalForce elementInternalForce = computeH8InternalForce(elementCoordinates ,elementDisplacements ,material);

            const H8TangentStiffness elementTangentStiffness = computeH8TangentStiffness(elementCoordinates ,elementDisplacements ,material);

            for(std::size_t j = 0 ;j < kH8DofCount ;++j) {
                result.internalForce[position[j]] += elementInternalForce[j];
            }

            for(std::size_t j = 0;j < kH8DofCount ;++j) {
                for(std::size_t k = 0;k < kH8DofCount ;++k) {
                    result.tangentStiffness[position[j]][position[k]] += elementTangentStiffness[j][k];
                }
            } 
        }
        return result;
    }

    [[nodiscard]]
    GlobalInternalForce_3D assembleInternalForce_3D (
        const H8Mesh& mesh,
        const GlobalDisplacement_3D& displacement,
        const MaterialParameters& material
    ){
        validateH8Mesh(mesh, displacement);

        const size_t N = mesh.nodes.size();

        GlobalInternalForce_3D internalForce{};

        internalForce.resize(3 * N, 0.0);


        for(std::size_t i = 0;i < mesh.elements.size() ;++i) {
            H8NodalValues elementCoordinates{};
            H8NodalValues elementDisplacements{};

            H8ElementConnectivity element = mesh.elements[i];

            std::array<std::size_t, kH8DofCount> position{};

            for(std::size_t j = 0 ;j < kH8NodeCount ;++j) {
                const size_t id = element[j];

                elementCoordinates[j] = mesh.nodes[id];

                elementDisplacements[j][0] = displacement[3 * id];
                elementDisplacements[j][1] = displacement[3 * id + 1];
                elementDisplacements[j][2] = displacement[3 * id + 2];

                position[3 * j] = 3 * id;
                position[3 * j + 1] = 3 * id + 1;
                position[3 * j + 2] = 3 * id + 2;
            }

            const H8InternalForce elementInternalForce = computeH8InternalForce(elementCoordinates ,elementDisplacements ,material);


            for(std::size_t j = 0 ;j < kH8DofCount ;++j) {
                internalForce[position[j]] += elementInternalForce[j];
            }
        }
        return internalForce;
    }

    [[nodiscard]]
    GlobalLumpedMass_3D assembleGlobalLumpedMass_3D (
        const H8Mesh& mesh,
        double referenceDensity
    ) {
        validateH8MeshGeometry(mesh);

        GlobalLumpedMass_3D mass{};

        const size_t N = mesh.nodes.size();

        mass.assign(3 * N,0.0);

        for(std::size_t i = 0;i < mesh.elements.size() ;++i) {
            H8NodalValues elementCoordinates{};
            H8ElementConnectivity element = mesh.elements[i];

            std::array<std::size_t, kH8DofCount> position{};

            for(std::size_t j = 0;j < kH8NodeCount;++j) {
                const size_t id = element[j];

                elementCoordinates[j] = mesh.nodes[id];

                position[3 * j] = 3 * id;
                position[3 * j + 1] = 3 * id + 1;
                position[3 * j + 2] = 3 * id + 2;
            }

            const H8LumpedMass h8Mass = computeH8LumpedMass(elementCoordinates,referenceDensity);

            for(std::size_t j = 0 ;j < kH8DofCount ;++j) {
                mass[position[j]] += h8Mass[j];
            }
        }
        return mass;
    }

    double assembleGlobalEnergy_3D (
        const H8Mesh& mesh,
        const GlobalDisplacement_3D& displacement,
        const MaterialParameters& material
    ){
        validateH8Mesh(mesh, displacement);

        const size_t N = mesh.elements.size();

        double energy = 0.0;

        for(std::size_t i = 0;i < N;++i) {
            H8NodalValues elementCoordinates{};
            H8NodalValues elementDisplacements{};

            H8ElementConnectivity element = mesh.elements[i];


            for(std::size_t j = 0 ;j < kH8NodeCount ;++j) {
                const size_t id = element[j];

                elementCoordinates[j] = mesh.nodes[id];

                elementDisplacements[j][0] = displacement[3 * id];
                elementDisplacements[j][1] = displacement[3 * id + 1];
                elementDisplacements[j][2] = displacement[3 * id + 2];
            }

            energy += computeH8StrainEnergy(elementCoordinates,elementDisplacements,material);
        }

        if (!std::isfinite(energy) || energy < 0.0) {
            throw std::invalid_argument(
                "Global strain energy must be finite and nonnegative."
            );
        }

        return energy;
    }

}