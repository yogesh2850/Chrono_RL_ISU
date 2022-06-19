// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Alessandro Tasora
// =============================================================================
//
// FEA (basic introduction)
//
// =============================================================================

//#include <chrono_pardisomkl/ChSolverPardisoMKL.h>
#include <solver/ChDirectSolverLS.h>
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/solver/ChSolverLS.h"
#include "chrono/fea/ChElementSpring.h"
#include "chrono/fea/ChElementTetra_4.h"
#include "chrono/fea/ChElementTetra_10.h"
#include "chrono/fea/ChElementHexa_8.h"
#include "chrono/fea/ChElementHexa_20.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/fea/ChElementBar.h"

#include "chrono/fea/ChContinuumThermal.h"

using namespace chrono;
using namespace fea;

// ====================================
// Test 1
// First example: SPRING ELEMENT
// ====================================

void loadMeshFromTetGenFile(std::shared_ptr<ChMesh> mesh,
                            const char* filename_node,
                            const char* filename_ele,
                            std::vector<std::shared_ptr<ChContinuumMaterial>> & my_material,
                            std::vector<int> & materialID) {
    int totnodes = 0;
    int nodes_offset = mesh->GetNnodes();
    int added_nodes = 0;

    // Load .node TetGen file
    {
        bool parse_header = true;
        bool parse_nodes = false;

        std::ifstream fin(filename_node);
        if (!fin.good())
            throw ChException("ERROR opening TetGen .node file: " + std::string(filename_node) + "\n");

        int nnodes = 0;
        int ndims = 0;
        int nattrs = 0;
        int nboundarymark = 0;

        std::string line;
        while (std::getline(fin, line)) {
            // trims white space from the beginning of the string
            line.erase(line.begin(),
                       std::find_if(line.begin(), line.end(), [](unsigned char c) { return !std::isspace(c); }));

            if (line[0] == '#')
                continue;  // skip comment
            if (line[0] == 0)
                continue;  // skip empty lines

            if (parse_header) {
                std::stringstream(line) >> nnodes >> ndims >> nattrs >> nboundarymark;
                if (ndims != 3)
                    throw ChException("ERROR in TetGen .node file. Only 3 dimensional nodes supported: \n" + line);
                if (nattrs != 0)
                    throw ChException("ERROR in TetGen .node file. Only nodes with 0 attrs supported: \n" + line);
                if (nboundarymark != 0)
                    throw ChException("ERROR in TetGen .node file. Only nodes with 0 markers supported: \n" + line);
                parse_header = false;
                parse_nodes = true;
                totnodes = nnodes;
                continue;
            }

            int idnode = 0;
            double x = -10e30;
            double y = -10e30;
            double z = -10e30;

            if (parse_nodes) {
                std::stringstream(line) >> idnode >> x >> y >> z;
                ++added_nodes;
                if (idnode <= 0 || idnode > nnodes)
                    throw ChException("ERROR in TetGen .node file. Node ID not in range: \n" + line + "\n");
                if (idnode != added_nodes)
                    throw ChException("ERROR in TetGen .node file. Nodes IDs must be sequential (1 2 3 ..): \n" + line +
                                      "\n");
                if (x == -10e30 || y == -10e30 || z == -10e30)
                    throw ChException("ERROR in TetGen .node file, in parsing x,y,z coordinates of node: \n" + line +
                                      "\n");

                ChVector<> node_position(x, y, z);

                if (true) {
                    auto mnode = chrono_types::make_shared<ChNodeFEAxyz>(node_position);
                    mesh->AddNode(mnode);
                } else
                    throw ChException("ERROR in TetGen generation. Material type not supported. \n");
            }

        }  // end while

    }  // end .node file

    // Load .ele TetGen file
    {
        bool parse_header = true;
        bool parse_tet = false;

        std::ifstream fin(filename_ele);
        if (!fin.good())
            throw ChException("ERROR opening TetGen .node file: " + std::string(filename_node) + "\n");

        int ntets, nnodespertet, nattrs = 0;

        std::string line;
        while (std::getline(fin, line)) {
            // trims white space from the beginning of the string
            line.erase(line.begin(),
                       std::find_if(line.begin(), line.end(), [](unsigned char c) { return !std::isspace(c); }));

            if (line[0] == '#')
                continue;  // skip comment
            if (line[0] == 0)
                continue;  // skip empty lines

            if (parse_header) {
                std::stringstream(line) >> ntets >> nnodespertet >> nattrs;
                if (nnodespertet != 4)
                    throw ChException("ERROR in TetGen .ele file. Only 4 -nodes per tes supported: \n" + line + "\n");
                if (nattrs != 0)
                    throw ChException("ERROR in TetGen .ele file. Only tets with 0 attrs supported: \n" + line + "\n");
                parse_header = false;
                parse_tet = true;
                continue;
            }

            int idtet = 0;
            int n1, n2, n3, n4;

            if (parse_tet) {
                std::stringstream(line) >> idtet >> n1 >> n2 >> n3 >> n4;
                if (idtet <= 0 || idtet > ntets)
                    throw ChException("ERROR in TetGen .node file. Tetrahedron ID not in range: \n" + line + "\n");
                if (n1 > totnodes)
                    throw ChException("ERROR in TetGen .node file, ID of 1st node is out of range: \n" + line + "\n");
                if (n2 > totnodes)
                    throw ChException("ERROR in TetGen .node file, ID of 2nd node is out of range: \n" + line + "\n");
                if (n3 > totnodes)
                    throw ChException("ERROR in TetGen .node file, ID of 3rd node is out of range: \n" + line + "\n");
                if (n4 > totnodes)
                    throw ChException("ERROR in TetGen .node file, ID of 4th node is out of range: \n" + line + "\n");

                int matID = 1;
                if(std::binary_search(materialID.begin(),materialID.end(),idtet - 1)){
                    matID = 0;
                }
                if (std::dynamic_pointer_cast<ChContinuumElastic>(my_material[matID])) {
                    auto mel = chrono_types::make_shared<ChElementTetra_4>();
                    mel->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(nodes_offset + n1 - 1)),
                                  std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(nodes_offset + n3 - 1)),
                                  std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(nodes_offset + n2 - 1)),
                                  std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(nodes_offset + n4 - 1)));
                    mel->SetMaterial(std::static_pointer_cast<ChContinuumElastic>(my_material[matID]));
                    mesh->AddElement(mel);
//                    std::cout << idtet << " " << my_material[matID]->Get_density() << "\n";
                }
                else
                    throw ChException("ERROR in TetGen generation. Material type not supported. \n");
            }

        }  // end while

    }  // end .ele file
}
void test_1() {
    GetLog() << "\n-------------------------------------------------\n";
    GetLog() << "TEST: spring element FEM  \n\n";

    // The physical system: it contains all physical objects.
    ChSystemSMC my_system;

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create some nodes. These are the classical point-like
    // nodes with x,y,z degrees of freedom, that can be used
    // for many types of FEM elements in space.
    // While creating them, also set X0 undeformed positions.
    auto mnodeA = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 0));
    auto mnodeB = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 1, 0));

    // For example, you can attach local 'point masses' (FE node masses are zero by default)
    mnodeA->SetMass(0.01);
    mnodeB->SetMass(0.01);

    // For example, set an applied force to a node:
    mnodeB->SetForce(ChVector<>(0, 5, 0));

    // Remember to add nodes and elements to the mesh!
    my_mesh->AddNode(mnodeA);
    my_mesh->AddNode(mnodeB);

    // Create some elements of 'spring-damper' type, each connecting
    // two 3D nodes:
    auto melementA = chrono_types::make_shared<ChElementSpring>();
    melementA->SetNodes(mnodeA, mnodeB);
    melementA->SetSpringK(100000);

    // Remember to add elements to the mesh!
    my_mesh->AddElement(melementA);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    // Create also a truss
    auto truss = chrono_types::make_shared<ChBody>();
    truss->SetBodyFixed(true);
    my_system.Add(truss);

    // Create a constraint between a node and the truss
    auto constraintA = chrono_types::make_shared<ChLinkPointFrame>();

    constraintA->Initialize(mnodeA,  // node to connect
                            truss);  // body to be connected to

    my_system.Add(constraintA);

    // Set no gravity
    // my_system.Set_G_acc(VNULL);

    // Perform a linear static analysis
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    my_system.SetSolver(solver);
    solver->SetMaxIterations(40);
    solver->SetTolerance(1e-10);
    solver->EnableDiagonalPreconditioner(true);
    solver->SetVerbose(true);

    my_system.DoStaticLinear();

    // Output result
    GetLog() << "poss after linear static analysis: \n";
    GetLog() << "  nodeA->pos \n" << mnodeA->GetPos();
    GetLog() << "  nodeB->pos \n" << mnodeB->GetPos();
    GetLog() << "Forces after linear static analysis: \n";
    GetLog() << "  constraintA.react \n" << constraintA->GetReactionOnBody();
}

//////////////////////////////////////////////////////////////////
// ============================================================ //
// Test 2													    //
// Second example: LINEAR TETRAHEDRAL ELEMENT				    //
// ============================================================ //
void test_2() {
    GetLog() << "\n-------------------------------------------------\n";
    GetLog() << "TEST: LINEAR tetrahedral element FEM  \n\n";

    // The physical system: it contains all physical objects.
    ChSystemSMC my_system;

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create a material, that must be assigned to each element,
    // and set its parameters
    auto mmaterial1 = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial1->Set_E(0.01e9);  // rubber 0.01e9, steel 200e9
    mmaterial1->Set_v(0.3);
    auto mmaterial2 = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial1->Set_E(200e9);  // rubber 0.01e9, steel 200e9
    mmaterial1->Set_v(0.3);

    // Create some nodes. These are the classical point-like
    // nodes with x,y,z degrees of freedom, that can be used
    // for many types of FEM elements in space.
    // While creating them, also set X0 undeformed positions.
    auto mnode1 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 0));
    auto mnode2 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 1));
    auto mnode3 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 1, 0));
    auto mnode4 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(1, 0, 0));

    auto mnode5 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(5 + 0, 5 + 0, 5 + 0));
    auto mnode6 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(5 + 0, 5 + 0, 5 + 1));
    auto mnode7 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(5 + 0, 5 + 1, 5 + 0));
    auto mnode8 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(5 + 1, 5 + 0, 5 + 0));

    // For example, set an applied force to a node:
    mnode3->SetForce(ChVector<>(0, 10000, 0));

    // Remember to add nodes and elements to the mesh!
    my_mesh->AddNode(mnode1);
    my_mesh->AddNode(mnode2);
    my_mesh->AddNode(mnode3);
    my_mesh->AddNode(mnode4);

    my_mesh->AddNode(mnode5);
    my_mesh->AddNode(mnode6);
    my_mesh->AddNode(mnode7);
    my_mesh->AddNode(mnode8);

    // Create the tetrahedron element, and assign
    // nodes and material
    auto melement1 = chrono_types::make_shared<ChElementTetra_4>();
    melement1->SetNodes(mnode1, mnode2, mnode3, mnode4);
    melement1->SetMaterial(mmaterial1);

    auto melement2 = chrono_types::make_shared<ChElementTetra_4>();
    melement2->SetNodes(mnode5, mnode6, mnode7, mnode8);
    melement2->SetMaterial(mmaterial2);

    // Remember to add elements to the mesh!
    my_mesh->AddElement(melement1);
    my_mesh->AddElement(melement2);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    // Create also a truss
    auto truss = chrono_types::make_shared<ChBody>();
    truss->SetBodyFixed(true);
    my_system.Add(truss);

    // Create a constraint between a node and the truss
    auto constraint1 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint2 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint3 = chrono_types::make_shared<ChLinkPointFrame>();

    auto constraint4 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint5 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint6 = chrono_types::make_shared<ChLinkPointFrame>();

    constraint1->Initialize(mnode1,  // node
                            truss);  // body to be connected to

    constraint2->Initialize(mnode2,  // node
                            truss);  // body to be connected to

    constraint3->Initialize(mnode4,  // node
                            truss);  // body to be connected to

    constraint4->Initialize(mnode5,  // node
                            truss);  // body to be connected to

    constraint5->Initialize(mnode6,  // node
                            truss);  // body to be connected to

    constraint6->Initialize(mnode8,  // node
                            truss);  // body to be connected to

    my_system.Add(constraint1);
    my_system.Add(constraint2);
    my_system.Add(constraint3);

    my_system.Add(constraint4);
    my_system.Add(constraint5);
    my_system.Add(constraint6);

    //    my_system.DumpSystemMatrices(true, true, true, true,"/home/maksbh/");
    // Set no gravity
    // my_system.Set_G_acc(VNULL);

    // Perform a linear static analysis
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    my_system.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-10);
    solver->EnableDiagonalPreconditioner(true);
    solver->SetVerbose(true);

    my_system.DoStaticLinear();

    // Output result
    GetLog() << "Resulting node positions:\n";
    GetLog() << mnode1->pos << "\n";
    GetLog() << mnode2->pos << "\n";
    GetLog() << mnode3->pos << "\n";
    GetLog() << mnode4->pos << "\n";

    GetLog() << mnode5->pos << "\n";
    GetLog() << mnode6->pos << "\n";
    GetLog() << mnode7->pos << "\n";
    GetLog() << mnode8->pos << "\n";

    GetLog() << "Resulting constraint reactions:\n";
    GetLog() << constraint1->GetReactionOnBody();
    GetLog() << constraint2->GetReactionOnBody();
    GetLog() << constraint3->GetReactionOnBody();
    GetLog() << constraint4->GetReactionOnBody();
    GetLog() << constraint5->GetReactionOnBody();
    GetLog() << constraint6->GetReactionOnBody();
}

//////////////////////////////////////////////////////////////////
// ============================================================ //
// Test 3													    //
// Second example: QUADRATIC TETRAHEDRAL ELEMENT				//
// ============================================================ //
void test_3() {
    GetLog() << "\n-------------------------------------------------\n";
    GetLog() << "TEST: QUADRATIC tetrahedral element FEM  \n\n";

    // The physical system: it contains all physical objects.
    ChSystemSMC my_system;

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create a material, that must be assigned to each element,
    // and set its parameters
    auto mmaterial = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial->Set_E(207e9);
    mmaterial->Set_v(0.3);

    // Create some nodes. These are the classical point-like
    // nodes with x,y,z degrees of freedom, that can be used
    // for many types of FEM elements in space.
    // While creating them, also set X0 undeformed positions.
    auto mnode1 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 0));
    auto mnode2 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0.001, 0, 0));
    auto mnode3 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0.001, 0));
    auto mnode4 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 0.001));
    auto mnode5 =
        chrono_types::make_shared<ChNodeFEAxyz>((mnode1->pos + mnode2->pos) * 0.5);  //  nodes at mid length of edges
    auto mnode6 = chrono_types::make_shared<ChNodeFEAxyz>((mnode2->pos + mnode3->pos) * 0.5);
    auto mnode7 = chrono_types::make_shared<ChNodeFEAxyz>((mnode3->pos + mnode1->pos) * 0.5);
    auto mnode8 = chrono_types::make_shared<ChNodeFEAxyz>((mnode1->pos + mnode4->pos) * 0.5);
    auto mnode9 = chrono_types::make_shared<ChNodeFEAxyz>((mnode4->pos + mnode2->pos) * 0.5);
    auto mnode10 = chrono_types::make_shared<ChNodeFEAxyz>((mnode3->pos + mnode4->pos) * 0.5);

    // For example, set an applied force to a node:
    mnode3->SetForce(ChVector<>(0, -1000, 0));

    // Remember to add nodes and elements to the mesh!
    my_mesh->AddNode(mnode1);
    my_mesh->AddNode(mnode2);
    my_mesh->AddNode(mnode3);
    my_mesh->AddNode(mnode4);
    my_mesh->AddNode(mnode5);
    my_mesh->AddNode(mnode6);
    my_mesh->AddNode(mnode7);
    my_mesh->AddNode(mnode8);
    my_mesh->AddNode(mnode9);
    my_mesh->AddNode(mnode10);

    // Create the tetrahedron element, and assign
    // it nodes and material
    auto melement1 = chrono_types::make_shared<ChElementTetra_10>();
    melement1->SetNodes(mnode1, mnode2, mnode3, mnode4, mnode5, mnode6, mnode7, mnode8, mnode9, mnode10);
    melement1->SetMaterial(mmaterial);

    // Remember to add elements to the mesh!
    my_mesh->AddElement(melement1);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    // Create also a truss
    auto truss = chrono_types::make_shared<ChBody>();
    my_system.Add(truss);
    truss->SetBodyFixed(true);

    // Create a constraint between a node and the truss
    auto constraint1 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint2 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint3 = chrono_types::make_shared<ChLinkPointFrame>();

    constraint1->Initialize(mnode1,  // node
                            truss);  // body to be connected to

    constraint2->Initialize(mnode2,  // node
                            truss);  // body to be connected to

    constraint3->Initialize(mnode4,  // node
                            truss);  // body to be connected to

    my_system.Add(constraint1);
    my_system.Add(constraint2);
    my_system.Add(constraint3);

    // Set no gravity
    // my_system.Set_G_acc(VNULL);

    // Perform a linear static analysis
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    my_system.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-12);
    solver->EnableDiagonalPreconditioner(true);
    solver->SetVerbose(true);

    my_system.DoStaticLinear();

    // Output result
    // GetLog()<<melement1.GetStiffnessMatrix()<<"\n";
    // GetLog()<<melement1.GetMatrB()<<"\n";
    GetLog() << mnode1->GetPos() << "\n";
    GetLog() << mnode2->GetPos() << "\n";
    GetLog() << mnode3->GetPos() << "\n";
    GetLog() << mnode4->GetPos() << "\n";
    GetLog() << "node3 displ: " << mnode3->GetPos() - mnode3->GetX0() << "\n";
}

//////////////////////////////////////////////////////////////////
// ============================================================ //
// Test 4													    //
// Second example: LINEAR HEXAHEDRAL ELEMENT					//
// ============================================================ //
void test_4() {
    GetLog() << "\n-------------------------------------------------\n";
    GetLog() << "TEST: LINEAR hexahedral element FEM  \n\n";

    // The physical system: it contains all physical objects.
    ChSystemSMC my_system;

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create a material, that must be assigned to each element,
    // and set its parameters
    auto mmaterial = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial->Set_E(207e6);
    mmaterial->Set_v(0.3);

    // Create some nodes. These are the classical point-like
    // nodes with x,y,z degrees of freedom, that can be used
    // for many types of FEM elements in space.
    // While creating them, also set X0 undeformed positions.
    double sx = 0.01;
    double sy = 0.10;
    double sz = 0.01;
    auto mnode1 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 0));
    auto mnode2 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, sz));
    auto mnode3 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(sx, 0, sz));
    auto mnode4 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(sx, 0, 0));
    auto mnode5 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, sy, 0));
    auto mnode6 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, sy, sz));
    auto mnode7 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(sx, sy, sz));
    auto mnode8 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(sx, sy, 0));

    // For example, set applied forces to nodes:
    mnode5->SetForce(ChVector<>(0, -1000, 0));
    mnode6->SetForce(ChVector<>(0, -1000, 0));
    mnode7->SetForce(ChVector<>(0, -1000, 0));
    mnode8->SetForce(ChVector<>(0, -1000, 0));

    // Remember to add nodes and elements to the mesh!
    my_mesh->AddNode(mnode1);
    my_mesh->AddNode(mnode2);
    my_mesh->AddNode(mnode3);
    my_mesh->AddNode(mnode4);
    my_mesh->AddNode(mnode5);
    my_mesh->AddNode(mnode6);
    my_mesh->AddNode(mnode7);
    my_mesh->AddNode(mnode8);

    // Create the tetrahedron element, and assign
    // it nodes and material
    auto melement1 = chrono_types::make_shared<ChElementHexa_8>();
    melement1->SetNodes(mnode1, mnode2, mnode3, mnode4, mnode5, mnode6, mnode7, mnode8);
    melement1->SetMaterial(mmaterial);

    // Remember to add elements to the mesh!
    my_mesh->AddElement(melement1);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    // Create also a truss
    auto truss = chrono_types::make_shared<ChBody>();
    my_system.Add(truss);
    truss->SetBodyFixed(true);

    // Create a constraint between a node and the truss
    auto constraint1 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint2 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint3 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint4 = chrono_types::make_shared<ChLinkPointFrame>();

    constraint1->Initialize(mnode1,  // node
                            truss);  // body to be connected to

    constraint2->Initialize(mnode2,  // node
                            truss);  // body to be connected to

    constraint3->Initialize(mnode3,  // node
                            truss);  // body to be connected to

    constraint4->Initialize(mnode4,  // node
                            truss);  // body to be connected to

    my_system.Add(constraint1);
    my_system.Add(constraint2);
    my_system.Add(constraint3);
    my_system.Add(constraint4);

    // Set no gravity
    // my_system.Set_G_acc(VNULL);

    // Perform a linear static analysis
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    my_system.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-12);
    solver->EnableDiagonalPreconditioner(true);
    solver->SetVerbose(true);

    my_system.DoStaticLinear();

    // Output result
    // GetLog()<<melement1.GetStiffnessMatrix()<<"\n";
    // GetLog()<<melement1.GetMatrB()<<"\n";
    GetLog() << mnode1->GetPos() << "\n";
    GetLog() << mnode2->GetPos() << "\n";
    GetLog() << mnode3->GetPos() << "\n";
    GetLog() << mnode4->GetPos() << "\n";
    GetLog() << "node5 displ: " << mnode5->GetPos() - mnode5->GetX0() << "\n";
    GetLog() << "node6 displ: " << mnode6->GetPos() - mnode6->GetX0() << "\n";
    GetLog() << "node7 displ: " << mnode7->GetPos() - mnode7->GetX0() << "\n";
    GetLog() << "node8 displ: " << mnode8->GetPos() - mnode8->GetX0() << "\n";
}

//////////////////////////////////////////////////////////////////
// ============================================================ //
// Test 5													    //
// Second example: QUADRATIC HEXAHEDRAL ELEMENT					//
// ============================================================ //
void test_5() {
    GetLog() << "\n-------------------------------------------------\n";
    GetLog() << "TEST: QUADRATIC hexahedral element FEM  \n\n";

    // The physical system: it contains all physical objects.
    ChSystemSMC my_system;

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create a material, that must be assigned to each element,
    // and set its parameters
    auto mmaterial = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial->Set_E(207e6);
    mmaterial->Set_v(0.3);

    // Create some nodes. These are the classical point-like
    // nodes with x,y,z degrees of freedom, that can be used
    // for many types of FEM elements in space.
    // While creating them, also set X0 undeformed positions.
    double sx = 0.01;
    double sy = 0.1;
    double sz = 0.01;
    auto mnode1 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 0));
    auto mnode2 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, sz));
    auto mnode3 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(sx, 0, sz));
    auto mnode4 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(sx, 0, 0));
    auto mnode5 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, sy, 0));
    auto mnode6 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, sy, sz));
    auto mnode7 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(sx, sy, sz));
    auto mnode8 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(sx, sy, 0));
    auto mnode9 = chrono_types::make_shared<ChNodeFEAxyz>((mnode1->pos + mnode2->pos) * 0.5);  // in between front face
    auto mnode10 = chrono_types::make_shared<ChNodeFEAxyz>((mnode2->pos + mnode3->pos) * 0.5);
    auto mnode11 = chrono_types::make_shared<ChNodeFEAxyz>((mnode3->pos + mnode4->pos) * 0.5);
    auto mnode12 = chrono_types::make_shared<ChNodeFEAxyz>((mnode1->pos + mnode4->pos) * 0.5);
    auto mnode13 = chrono_types::make_shared<ChNodeFEAxyz>((mnode5->pos + mnode6->pos) * 0.5);  // in between back face
    auto mnode14 = chrono_types::make_shared<ChNodeFEAxyz>((mnode6->pos + mnode7->pos) * 0.5);
    auto mnode15 = chrono_types::make_shared<ChNodeFEAxyz>((mnode7->pos + mnode8->pos) * 0.5);
    auto mnode16 = chrono_types::make_shared<ChNodeFEAxyz>((mnode8->pos + mnode5->pos) * 0.5);
    auto mnode17 = chrono_types::make_shared<ChNodeFEAxyz>((mnode2->pos + mnode6->pos) * 0.5);  // in between side edges
    auto mnode18 = chrono_types::make_shared<ChNodeFEAxyz>((mnode3->pos + mnode7->pos) * 0.5);
    auto mnode19 = chrono_types::make_shared<ChNodeFEAxyz>((mnode4->pos + mnode8->pos) * 0.5);
    auto mnode20 = chrono_types::make_shared<ChNodeFEAxyz>((mnode1->pos + mnode5->pos) * 0.5);

    // For example, set applied forces to nodes:
    mnode5->SetForce(ChVector<>(0, -500, 0));
    mnode6->SetForce(ChVector<>(0, -500, 0));
    mnode7->SetForce(ChVector<>(0, -500, 0));
    mnode8->SetForce(ChVector<>(0, -500, 0));
    mnode13->SetForce(ChVector<>(0, -500, 0));
    mnode14->SetForce(ChVector<>(0, -500, 0));
    mnode15->SetForce(ChVector<>(0, -500, 0));
    mnode16->SetForce(ChVector<>(0, -500, 0));

    // Remember to add nodes and elements to the mesh!
    my_mesh->AddNode(mnode1);
    my_mesh->AddNode(mnode2);
    my_mesh->AddNode(mnode3);
    my_mesh->AddNode(mnode4);
    my_mesh->AddNode(mnode5);
    my_mesh->AddNode(mnode6);
    my_mesh->AddNode(mnode7);
    my_mesh->AddNode(mnode8);
    my_mesh->AddNode(mnode9);
    my_mesh->AddNode(mnode10);
    my_mesh->AddNode(mnode11);
    my_mesh->AddNode(mnode12);
    my_mesh->AddNode(mnode13);
    my_mesh->AddNode(mnode14);
    my_mesh->AddNode(mnode15);
    my_mesh->AddNode(mnode16);
    my_mesh->AddNode(mnode17);
    my_mesh->AddNode(mnode18);
    my_mesh->AddNode(mnode19);
    my_mesh->AddNode(mnode20);

    // Create the tetrahedron element, and assign
    // its nodes and material
    auto melement1 = chrono_types::make_shared<ChElementHexa_20>();
    melement1->SetNodes(mnode1, mnode2, mnode3, mnode4, mnode5, mnode6, mnode7, mnode8, mnode9, mnode10, mnode11,
                        mnode12, mnode13, mnode14, mnode15, mnode16, mnode17, mnode18, mnode19, mnode20);
    melement1->SetMaterial(mmaterial);

    // Use this statement to use the reduced integration
    // Default number of gauss point: 27. Reduced integration -> 8 Gp.
    melement1->SetReducedIntegrationRule();

    // Remember to add elements to the mesh!
    my_mesh->AddElement(melement1);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    // Create also a truss
    auto truss = chrono_types::make_shared<ChBody>();
    my_system.Add(truss);
    truss->SetBodyFixed(true);

    // Create a constraint between a node and the truss
    auto constraint1 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint2 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint3 = chrono_types::make_shared<ChLinkPointFrame>();
    auto constraint4 = chrono_types::make_shared<ChLinkPointFrame>();

    constraint1->Initialize(mnode1,  // node
                            truss);  // body to be connected to

    constraint2->Initialize(mnode2,  // node
                            truss);  // body to be connected to

    constraint3->Initialize(mnode3,  // node
                            truss);  // body to be connected to

    constraint4->Initialize(mnode4,  // node
                            truss);  // body to be connected to

    my_system.Add(constraint1);
    my_system.Add(constraint2);
    my_system.Add(constraint3);
    my_system.Add(constraint4);

    // Set no gravity
    // my_system.Set_G_acc(VNULL);

    // Perform a linear static analysis
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    my_system.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-12);
    solver->EnableDiagonalPreconditioner(true);
    solver->SetVerbose(true);

    my_system.DoStaticLinear();

    // Output some results
    GetLog() << "node5 displ: " << mnode5->GetPos() - mnode5->GetX0() << "\n";
    GetLog() << "node6 displ: " << mnode6->GetPos() - mnode6->GetX0() << "\n";
    GetLog() << "node7 displ: " << mnode7->GetPos() - mnode7->GetX0() << "\n";
    GetLog() << "node8 displ: " << mnode8->GetPos() - mnode8->GetX0() << "\n";
    GetLog() << "Element volume" << melement1->GetVolume() << "\n";
}

// Do some tests in a single run, inside the main() function.
// Results will be simply text-formatted outputs in the console..

void dumpData(const ChMesh& chmesh, const std::string& fileName) {
    static constexpr int npes = 4;
    int numCells = chmesh.GetElements().size();
    int numVertices = numCells * npes;
    FILE* fp = fopen(fileName.c_str(), "w+");
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
    fprintf(fp, "<UnstructuredGrid >\n");
    fprintf(fp, "<Piece NumberOfPoints=\" %d \" NumberOfCells=\" %d \" >\n", numVertices, numCells);
    { /** Points data **/
        float* coords = new float[numVertices * 3];
        /** Parview expects point to be written in 3D**/
        fprintf(fp, "<Points>\n");
        fprintf(fp, "<DataArray type=\"Float32\" NumberOfComponents=\" 3\" format=\"ascii\">\n");
        for (int i = 0; i < numCells; i++) {
            for (int numELe = 0; numELe < npes; numELe++) {
                auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(chmesh.GetElements()[i]->GetNodeN(numELe));
                fprintf(fp, "%f %f %f\n", mnode->GetPos().x(), mnode->GetPos().y(), mnode->GetPos().z());
                //                std::cout << "In\n";
            }
        }
        delete[] coords;
        fprintf(fp, "\n");
        fprintf(fp, "</DataArray>\n");
        fprintf(fp, "</Points>\n");
    }
    fprintf(fp, "<Cells>\n");
    {
        /** Connectivity **/

        fprintf(fp, "<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n");
        for (unsigned int ele = 0; ele < numCells; ele++) {
            fprintf(fp, "%d %d %d %d\n", ele * npes, ele * npes + 1, ele * npes + 2, ele * npes + 3);
        }
        fprintf(fp, "</DataArray>\n");
    }
    {
        /** Offsets **/
        fprintf(fp, "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n");
        for (unsigned int ele = 0; ele < numCells; ele++) {
            fprintf(fp, "%d\n", 4 * (ele + 1));
        }
        fprintf(fp, "</DataArray>\n");
    }
    {
        /** Write cellTypes **/
        fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
        for (unsigned int ele = 0; ele < numCells; ele++) {
            fprintf(fp, "%d ", 10);
        }
        fprintf(fp, "</DataArray>\n");
    }
    fprintf(fp, "</Cells>\n");

    //    {
    //        fprintf(fp, "<PointData>\n");
    //        fprintf(fp, "<DataArray type=\"Float64\" Name=\" %s \" format=\"ascii\">\n", "Temperature");
    //        for(int i = 0; i < numCells; i++){
    //            for(int numELe = 0; numELe < npes; numELe++) {
    //                auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyzP>(chmesh.GetElements()[i]->GetNodeN(numELe));
    //                fprintf(fp,"%f\n",mnode->GetP());
    //            }
    //        }
    //        fprintf(fp, "</DataArray>\n");
    //        fprintf(fp, "</PointData>\n");
    //    }
    fprintf(fp, "</Piece>\n");
    fprintf(fp, "</UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}

void solveElasticity() {
    GetLog() << "\n-------------------------------------------------\n";
    GetLog() << "TEST: tetrahedron FEM dynamics, implicit integration \n\n";

    // The physical system: it contains all physical objects.
    ChSystemSMC my_system;

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create a material, that must be assigned to each element,
    // and set its parameters
    auto mmaterial = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial->Set_E(200e9);  // rubber 0.01e9, steel 200e9
    mmaterial->Set_v(0.3);
    mmaterial->Set_RayleighDampingK(0.01);
    mmaterial->Set_density(1000);

    int numBoundary;
    std::ifstream fin("/media/maksbh/1e2e184f-abc8-40c0-b27c-a22673998855/run/chrono/nodeID.txt", std::ios::in);
    fin >> numBoundary;
    std::vector<int> boundaryIds(numBoundary);
    for (auto& bnd : boundaryIds) {
        fin >> bnd;
    }
    fin.close();

    ChMeshFileLoader::FromTetGenFile(
        my_mesh, "/media/maksbh/1e2e184f-abc8-40c0-b27c-a22673998855/run/chrono/maizeA.node",
        "/media/maksbh/1e2e184f-abc8-40c0-b27c-a22673998855/run/chrono/maizeA.ele", mmaterial);

    for (const auto& bnd : boundaryIds) {
        auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(bnd - 1));
        if (mnode->GetPos().z() > 10.0) {
            mnode->SetForce(ChVector<>(0, 0, 0));
        }
    }

    my_system.Add(my_mesh);
    auto truss = chrono_types::make_shared<ChBody>();
    truss->SetBodyFixed(true);
    my_system.Add(truss);

    // Create constraints between nodes and truss
    // (for example, fix to ground all nodes which are near y=0
    for (unsigned int inode = 0; inode < my_mesh->GetNnodes(); ++inode) {
        if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(inode))) {
            if (mnode->GetPos().z() < 2.0) {
                auto constraint = chrono_types::make_shared<ChLinkPointFrame>();
                constraint->Initialize(mnode, truss);
                my_system.Add(constraint);
            }
        }
    }

    // Create some nodes. These are the classical point-like
    // nodes with x,y,z degrees of freedom, that can be used
    // for many types of FEM elements in space.
    //    auto mnode1 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 0));
    //    auto mnode2 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0, 1));
    //    auto mnode3 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(0, 1, 0));
    //    auto mnode4 = chrono_types::make_shared<ChNodeFEAxyz>(ChVector<>(1, 0, 0));
    //
    //    // For example, set a point-like mass at a node:
    //    mnode1->SetMass(200);
    //    mnode2->SetMass(200);
    //    mnode3->SetMass(200);
    //    mnode4->SetMass(200);
    //    // For example, set an initial displacement to a node:
    //    mnode3->SetPos(mnode3->GetX0() + ChVector<>(0, 0.01, 0));
    //
    //    // Remember to add nodes and elements to the mesh!
    //    my_mesh->AddNode(mnode1);
    //    my_mesh->AddNode(mnode2);
    //    my_mesh->AddNode(mnode3);
    //    my_mesh->AddNode(mnode4);

    // Create the tetrahedron element, and assign
    //    // nodes and material
    //    auto melement1 = chrono_types::make_shared<ChElementTetra_4>();
    //    melement1->SetNodes(mnode1, mnode2, mnode3, mnode4);
    //    melement1->SetMaterial(mmaterial);

    // Remember to add elements to the mesh!
    //    my_mesh->AddElement(melement1);

    // Remember to add the mesh to the system!
    //    my_system.Add(my_mesh);
    //
    //    // Create also a truss
    //    auto truss = chrono_types::make_shared<ChBody>();
    //    truss->SetBodyFixed(true);
    //    my_system.Add(truss);
    //
    //    for (unsigned int inode = 0; inode < my_mesh->GetNnodes(); ++inode) {
    //        auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(inode));
    //        if(mnode->GetPos().z() < 0.1) {
    //            auto constraint1 = chrono_types::make_shared<ChLinkPointFrame>();
    //            constraint1->Initialize(mnode,  // node
    //                                    truss);  // body to be connected to
    //            my_system.Add(constraint1);
    //        }
    //    }
    // Create a constraint between a node and the truss
    //    auto constraint1 = chrono_types::make_shared<ChLinkPointFrame>();
    //    auto constraint2 = chrono_types::make_shared<ChLinkPointFrame>();
    //    auto constraint3 = chrono_types::make_shared<ChLinkPointFrame>();
    //
    //    constraint1->Initialize(mnode1,  // node
    //                            truss);  // body to be connected to
    //
    //    constraint2->Initialize(mnode2,  // node
    //                            truss);  // body to be connected to
    //
    //    constraint3->Initialize(mnode4,  // node
    //                            truss);  // body to be connected to

    //    my_system.Add(constraint1);
    //    my_system.Add(constraint2);
    //    my_system.Add(constraint3);

    // Perform a dynamic time integration:

    auto solver = chrono_types::make_shared<ChSolverSparseLU>();
    solver->UseSparsityPatternLearner(false);
    solver->LockSparsityPattern(false);
    solver->SetVerbose(true);
    my_system.SetSolver(solver);

    //    auto solver = chrono_types::make_shared<ChSolverGMRES>();
    //   my_system.SetSolver(solver);
    //   solver->SetMaxIterations(1000);
    //   solver->SetTolerance(1e-5);
    //   solver->EnableDiagonalPreconditioner(true);
    // solver->SetVerbose(true);
    //  solver->EnableWarmStart(true);

    my_system.SetSolverForceTolerance(1e-10);

    my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);

    int counter = 0;
    double timestep = 0.1;
    while (my_system.GetChTime() < 4) {
        my_system.DoStepDynamics(timestep);
        dumpData(*my_mesh, "elastic" + std::to_string(counter) + ".vtu");
        counter++;
        if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(my_mesh->GetNnodes() - 1))) {
            GetLog() << 0 << " " << mnode->GetPos().x() << " " << mnode->GetPos().y() << " " << mnode->GetPos().z()
                     << "\n";
        }
    }
}

void solveVariableElasticity(const std::string & boundaryFileName, const std::string & stemCells,
                             const std::string & nodeList, const std::string & connectivityList) {
    // The physical system: it contains all physical objects.
    ChSystemSMC my_system;
    ChVector<> acc = {0,0,-9.8};
    my_system.Set_G_acc(acc);
    // Material properties for material 1 (Stem)
    double E1 = 16e9;// rubber 0.01e9, steel 200e9
    double v1 = 0.3;
    double dampingCoeff1 = 0.1;
    double density1 = 400;

    // Material properties for material 2 (Leaf)
    double E2 = 0.01e9;// rubber 0.01e9, steel 200e9
    double v2 = 0.3;
    double dampingCoeff2 = 0.1; //make it higher (*10)
    double density2 = 10;

    double timestep = 1.0;
    double finalTime = 50;

    bool writeVTU = true;
    bool writeTXT = true;
    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

//    my_mesh->SetAutomaticGravity(true);

    auto mmaterial1 = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial1->Set_E(E1);
    mmaterial1->Set_v(v1);
    mmaterial1->Set_RayleighDampingK(dampingCoeff1);
    mmaterial1->Set_density(density1);

    auto mmaterial2 = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial2->Set_E(E2);  // rubber 0.01e9, steel 200e9
    mmaterial2->Set_v(v2);
    mmaterial2->Set_RayleighDampingK(dampingCoeff2);
    mmaterial2->Set_density(density2);

    std::vector<std::shared_ptr<ChContinuumMaterial>> material = {mmaterial1,mmaterial2};
    int numBoundary;
    struct BoundaryConstruct{
        int boundaryID;
        bool isDirichlet;
        double forceValue[3];
    };
    std::vector<BoundaryConstruct> boundaryIds;
    {
        std::ifstream fin(boundaryFileName, std::ios::in);
        fin >> numBoundary;
        boundaryIds.resize(numBoundary);
        for (auto& bnd : boundaryIds) {
            fin >> bnd.boundaryID;
            fin >> bnd.isDirichlet;
            if(not bnd.isDirichlet){
                fin >> bnd.forceValue[0];
                fin >> bnd.forceValue[1];
                fin >> bnd.forceValue[2];
            }

        }
        fin.close();
    }

    std::vector<int> materialID;
    {
        int numCells;
        std::ifstream fin(stemCells, std::ios::in);
        fin>> numCells;
        materialID.resize(numCells);
        for (auto& cell : materialID) {
            fin >> cell;
        }
        fin.close();
    }






    loadMeshFromTetGenFile(
        my_mesh, nodeList.c_str(),
        connectivityList.c_str(), material,materialID);

    my_system.Add(my_mesh);
    auto truss = chrono_types::make_shared<ChBody>();
    truss->SetBodyFixed(true);
    my_system.Add(truss);

    for (const auto& bnd : boundaryIds) {
        auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(bnd.boundaryID - 1));
        if(not(bnd.isDirichlet)) {
            mnode->SetForce(ChVector<>(bnd.forceValue[0], bnd.forceValue[1], bnd.forceValue[2]));
            if(bnd.forceValue[2]!=0)
            std::cout << mnode->GetPos().x() << " " << mnode->GetPos().y() << " " << mnode->GetPos().z() << "\n";
        }
//        else{
//            auto constraint = chrono_types::make_shared<ChLinkPointFrame>();
//            constraint->Initialize(mnode, truss);
//            my_system.Add(constraint);
//        }
    }




    for (unsigned int inode = 0; inode < my_mesh->GetNnodes(); ++inode) {

        if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(inode))) {
            if (mnode->GetPos().z() < 0) {
                auto constraint = chrono_types::make_shared<ChLinkPointFrame>();
                constraint->Initialize(mnode, truss);
                my_system.Add(constraint);
            }
        }
    }
    auto solver = chrono_types::make_shared<ChSolverSparseLU>();
    solver->UseSparsityPatternLearner(false);
    solver->LockSparsityPattern(false);
    solver->SetVerbose(false);
    my_system.SetSolver(solver);

    my_system.SetSolverForceTolerance(1e-10);

    my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);

    int counter = 0;

    while (my_system.GetChTime() < finalTime) {
        my_system.DoStepDynamics(timestep);
        if(writeVTU) {
            dumpData(*my_mesh, "elastic" + std::to_string(counter) + ".vtu");
        }
        counter++;
//        if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(my_mesh->GetNnodes() - 1))) {
//            //            if (mnode->GetPos().x() < 0.01) {
//            GetLog() << 0 << " " << mnode->GetPos().x() << " " << mnode->GetPos().y() << " " << mnode->GetPos().z()
//                     << "\n";
//            //            }
//        }
        //        GetLog() << " t =" << my_system.GetChTime() << "  mnode3 pos.y()=" << mnode3->GetPos().y() << "  \n";
    }
    if(writeTXT){
        std::ofstream fout("Position.txt");
        for (unsigned int inode = 0; inode < my_mesh->GetNnodes(); ++inode) {
            if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(inode))) {
                fout << mnode->GetPos().x() << " " <<mnode->GetPos().y() << " " << mnode->GetPos().z()<< "\n";
            }
        }
        fout.close();
    }


}

void solveHeat() {
    ChSystemSMC thermal;
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    auto mmaterial = chrono_types::make_shared<ChContinuumThermal>();
    mmaterial->SetMassSpecificHeatCapacity(2);
    mmaterial->SetThermalConductivityK(200);

    //    ChMeshFileLoader::FromTetGenFile(my_mesh,
    //    "/media/maksbh/1e2e184f-abc8-40c0-b27c-a22673998855/packages/chrono/data/fea/beam.node",
    //                                     "/media/maksbh/1e2e184f-abc8-40c0-b27c-a22673998855/packages/chrono/data/fea/beam.ele",mmaterial);
    ChMeshFileLoader::FromTetGenFile(
        my_mesh, "/media/maksbh/1e2e184f-abc8-40c0-b27c-a22673998855/run/chrono/maizeA.node",
        "/media/maksbh/1e2e184f-abc8-40c0-b27c-a22673998855/run/chrono/maizeA.ele", mmaterial);
    std::ifstream fin("/media/maksbh/1e2e184f-abc8-40c0-b27c-a22673998855/run/chrono/nodeID.txt", std::ios::in);

    int numBoundary;
    fin >> numBoundary;
    std::vector<int> boundaryIds(numBoundary);
    for (auto& bnd : boundaryIds) {
        fin >> bnd;
    }
    fin.close();
    auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyzP>(my_mesh->GetElements()[0]->GetNodeN(0));
    std::cout << mnode->GetPos().z() << "\n";
    //    std::fill(boundaryIds.begin(),boundaryIds.end(),1);

    for (const auto& bnd : boundaryIds) {
        auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyzP>(my_mesh->GetNode(bnd - 1));
        mnode->SetFixed(true);
        mnode->SetP(fabs(mnode->GetPos().z()));  // field: temperature [K]
        //        if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyzP>(my_mesh->GetNode(inode))) {
        //            if (mnode->GetPos().z() < 10) {
        //                mnode->SetFixed(true);
        //                mnode->SetP(10);  // field: temperature [K]
        //            }
        //            if (mnode->GetPos().z() > 70) {
        //                mnode->SetFixed(true);
        //                mnode->SetP(1);  // field: temperature [K]
        //            }
        //        }
    }
    thermal.Add(my_mesh);

    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    thermal.SetSolver(solver);
    solver->SetMaxIterations(500);
    solver->SetTolerance(1e-11);
    solver->EnableDiagonalPreconditioner(true);
    solver->EnableWarmStart(true);  // IMPORTANT for convergence when using EULER_IMPLICIT_LINEARIZED
    solver->SetVerbose(true);
    thermal.DoStaticLinear();
    for (unsigned int inode = 0; inode < my_mesh->GetNnodes(); ++inode) {
        auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyzP>(my_mesh->GetNode(inode));
        GetLog() << inode << " has T=" << mnode->GetP() << "\n";
    }

    dumpData(*my_mesh, "data.vtu");
}

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    GetLog() << " Example: the FEM technology for finite elements \n\n\n";

    if(argc < 5){
        std::cout << "Usage: BoundaryFileName stemCellsList nodeList ConnectivityList \n";
        exit(EXIT_FAILURE);
    }
    const std::string boundaryFileName = argv[1];
    const std::string stemCellsList = argv[2];
    const std::string nodeList = argv[3];
    const std::string ConnectivityList = argv[4];


    //    //test_1(); //// NOT WORKING
//    test_2();
    //    test_3();
    //    test_4();
    //    test_5();
    //    solveElasticity();
        solveVariableElasticity(boundaryFileName,stemCellsList,nodeList,ConnectivityList);
    //    solveHeat();

    return 0;
}
