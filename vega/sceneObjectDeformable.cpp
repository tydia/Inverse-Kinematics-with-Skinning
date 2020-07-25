/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "sceneObject" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Daniel Schroeder                         *
 * http://www.jernejbarbic.com/vega                                      *
 *                                                                       *
 * Research: Jernej Barbic, Hongyi Xu, Yijing Li,                        *
 *           Danyong Zhao, Bohan Wang,                                   *
 *           Fun Shing Sin, Daniel Schroeder,                            *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC,                *
 *          Sloan Foundation, Okawa Foundation,                          *
 *          USC Annenberg Foundation                                     *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sceneObjectDeformable.h"

SceneObjectDeformable::SceneObjectDeformable(const char * filenameOBJ):
   SceneObjectWithRestPosition(filenameOBJ) 
{
}

SceneObjectDeformable::SceneObjectDeformable(ObjMesh * objMesh, bool deepCopy):
   SceneObjectWithRestPosition(objMesh, deepCopy) 
{
}

SceneObjectDeformable::~SceneObjectDeformable()
{
}

void SceneObjectDeformable::ResetDeformationToRest()
{
  for(int i = 0; i < n; i++)
    mesh->setPosition(i, Vec3d(restPosition[3 * i + 0], restPosition[3 * i + 1], restPosition[3 * i + 2]));
}

void SceneObjectDeformable::AddVertexDeformations(const double * u)
{
  for(int i = 0; i < n; i++)
  {
    mesh->setPosition(i, mesh->getPosition(i) + Vec3d(u[3 * i + 0], u[3 * i + 1], u[3 * i + 2]));
  }
}

void SceneObjectDeformable::SetVertexDeformations(const double * u)
{
  for(int i = 0; i < n; i++)
  {
    mesh->setPosition(i, Vec3d(restPosition[3 * i + 0] + u[3 * i + 0], restPosition[3 * i + 1] + u[3 * i + 1], restPosition[3 * i + 2] + u[3 * i + 2]));
  }
}

void SceneObjectDeformable::SetVertexDeformations(const float * u)
{
  // set the deformations
  for(int i = 0; i < n; i++)
  {
    mesh->setPosition(i, Vec3d(restPosition[3 * i + 0] + u[3 * i + 0], restPosition[3 * i + 1] + u[3 * i + 1], restPosition[3 * i + 2] + u[3 * i + 2]));
  }
}

void SceneObjectDeformable::SetLighting(Lighting * lighting)
{
  lighting->LightScene();
}

