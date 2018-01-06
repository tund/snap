#ifndef N2VVEL_H
#define N2VVEL_H

#include "stdafx.h"

#include "Snap.h"
#include "velesrandomwalk.h"
#include "word2vec.h"

///Calculates node2vec feature representation for nodes and writes them into EmbeddinsHV
void Veles(PWNet& InNet, double& ParamP, double& ParamQ, int& Dimensions,
 int& WalkLen, int& NumWalks, int& WinSize, int& Iter, bool& Verbose,
 TIntFltVH& EmbeddingsHV, bool& OutWalks, TVVec<TInt, int64>& WalksVV); 

#endif //N2VVEL_H