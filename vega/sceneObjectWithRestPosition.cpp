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

#include <string.h>
#include <stdlib.h>
#include "sceneObjectWithRestPosition.h"

SceneObjectWithRestPosition::SceneObjectWithRestPosition(const char * filename): SceneObject(filename) 
{
  Construct();
}

SceneObjectWithRestPosition::SceneObjectWithRestPosition(ObjMesh * objMesh, bool deepCopy): SceneObject(objMesh, deepCopy) 
{
  Construct();
}

void SceneObjectWithRestPosition::Construct()
{
  restPosition = (double*) malloc (sizeof(double) * 3 * n);
  for(int i = 0; i < n; i++)
  {
    Vec3d pos = mesh->getPosition(i);
    restPosition[3 * i + 0] = pos[0];
    restPosition[3 * i + 1] = pos[1];
    restPosition[3 * i + 2] = pos[2];
  }
}

SceneObjectWithRestPosition::~SceneObjectWithRestPosition()
{
  free(restPosition);
}

void SceneObjectWithRestPosition::GetVertexRestPositions(double * buffer)
{
  memcpy(buffer, restPosition, 3 * n * sizeof(double));
}

void SceneObjectWithRestPosition::GetVertexRestPositions(float * buffer)
{
  for(int i=0; i<3*n; i++)
    buffer[i] = (float) (restPosition[i]);
}
