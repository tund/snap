#include "stdafx.h"
#include "Snap.h"
#include "biasedrandomwalk.h"

//Preprocess alias sampling method
void GetNodeAlias(TFltV& PTblV, TIntVFltVPr& NTTable) {
  int64 N = PTblV.Len();
  TIntV& KTbl = NTTable.Val1;
  TFltV& UTbl = NTTable.Val2;
  for (int64 i = 0; i < N; i++) {
    KTbl[i]=0;
    UTbl[i]=0;
  }
  TIntV UnderV;
  TIntV OverV;
  for (int64 i = 0; i < N; i++) {
    UTbl[i] = PTblV[i]*N;
    if (UTbl[i] < 1) {
      UnderV.Add(i);
    } else {
      OverV.Add(i);
    }
  }
  while (UnderV.Len() > 0 && OverV.Len() > 0) {
    int64 Small = UnderV.Last();
    int64 Large = OverV.Last();
    UnderV.DelLast();
    OverV.DelLast();
    KTbl[Small] = Large;
    UTbl[Large] = UTbl[Large] + UTbl[Small] - 1;
    if (UTbl[Large] < 1) {
      UnderV.Add(Large);
    } else {
      OverV.Add(Large);
    }
  }
}

//Get random element using alias sampling method
int64 AliasDrawInt(TIntVFltVPr& NTTable, TRnd& Rnd) {
  int64 N = NTTable.GetVal1().Len();
  TInt X = static_cast<int64>(Rnd.GetUniDev()*N);
  double Y = Rnd.GetUniDev();
  return Y < NTTable.GetVal2()[X] ? X : NTTable.GetVal1()[X];
}

void PreprocessNode (PWNet& InNet, double& ParamP, double& ParamQ,
 TWNet::TNodeI NI, int64& NCnt, bool& Verbose) {
  if (Verbose && NCnt%100 == 0) {
    printf("\rPreprocessing progress: %.2lf%% ",(double)NCnt*100/(double)(InNet->GetNodes()));fflush(stdout);
  }
  //for node t
  THash <TInt, TBool> NbrH;                                    //Neighbors of t
  for (int64 i = 0; i < NI.GetOutDeg(); i++) {
    NbrH.AddKey(NI.GetNbrNId(i));
  } 
  for (int64 i = 0; i < NI.GetOutDeg(); i++) {
    TWNet::TNodeI CurrI = InNet->GetNI(NI.GetNbrNId(i));      //for each node v
    double Psum = 0;
    TFltV PTable;                              //Probability distribution table
    for (int64 j = 0; j < CurrI.GetOutDeg(); j++) {           //for each node x
      int64 FId = CurrI.GetNbrNId(j);
      TFlt Weight;
      if (!(InNet->GetEDat(CurrI.GetId(), FId, Weight))){ continue; }
      if (FId==NI.GetId()) {
        PTable.Add(Weight / ParamP);
        Psum += Weight / ParamP;
      } else if (NbrH.IsKey(FId)) {
        PTable.Add(Weight);
        Psum += Weight;
      } else {
        PTable.Add(Weight / ParamQ);
        Psum += Weight / ParamQ;
      }
    }
    //Normalizing table
    for (int64 j = 0; j < CurrI.GetOutDeg(); j++) {
      PTable[j] /= Psum;
    }
    GetNodeAlias(PTable,CurrI.GetDat().GetDat(NI.GetId()));
  }
  NCnt++;
}

//Preprocess transition probabilities for each path t->v->x
void PreprocessTransitionProbs(PWNet& InNet, double& ParamP, double& ParamQ, bool& Verbose) {
  for (TWNet::TNodeI NI = InNet->BegNI(); NI < InNet->EndNI(); NI++) {
    InNet->SetNDat(NI.GetId(),TIntIntVFltVPrH());
  }
  printf("PreprocessTransitionProbs 1\n");
  for (TWNet::TNodeI NI = InNet->BegNI(); NI < InNet->EndNI(); NI++) {
    printf("%d\n", NI.GetId());
    for (int64 i = 0; i < NI.GetOutDeg(); i++) {                    //allocating space in advance to avoid issues with multithreading
      TWNet::TNodeI CurrI = InNet->GetNI(NI.GetNbrNId(i));
      CurrI.GetDat().AddDat(NI.GetId(),TPair<TIntV,TFltV>(TIntV(CurrI.GetOutDeg()),TFltV(CurrI.GetOutDeg())));
    }
  }
  int64 NCnt = 0;
  TIntV NIds;
  for (TWNet::TNodeI NI = InNet->BegNI(); NI < InNet->EndNI(); NI++) {
    NIds.Add(NI.GetId());
  }
#pragma omp parallel for schedule(dynamic)
  for (int64 i = 0; i < NIds.Len(); i++) {
    PreprocessNode(InNet, ParamP, ParamQ, InNet->GetNI(NIds[i]), NCnt, Verbose);
  }
  if(Verbose){ printf("\n"); }
}

int64 PredictMemoryRequirements(PWNet& InNet) {
  int count_node_out_deg = 0;
  for (TWNet::TNodeI NI = InNet->BegNI(); NI < InNet->EndNI(); NI++) {
    if (NI.GetOutDeg() >= 500) {
      count_node_out_deg += 1;
    }
  }
  printf("Num nodes with OutDeg >= 500: %d\n", count_node_out_deg);

  // FILE *fout;
  // fout = fopen("emb/edges.tmp", "w");
  

  long long int MemNeeded = 0;
  printf("BegNI.GetId, EndNI.GetId: %d\t%d\n", InNet->BegNI().GetId(), InNet->EndNI().GetId());
  int count_node = 0;
  int count_edge = 0;
  int max_out_deg = 0;
  int min_out_deg = 1000000000;
  int64 sum_out_deg = 0;
  int64 count_200k = 0;
  printf("sizeof(TInt) + sizeof(TFlt) = %d", sizeof(TInt) + sizeof(TFlt));
  for (TWNet::TNodeI NI = InNet->BegNI(); NI < InNet->EndNI(); NI++) {
    count_node += 1;
    for (int64 i = 0; i < NI.GetOutDeg(); i++) {
      count_edge += 1;
      TWNet::TNodeI CurrI = InNet->GetNI(NI.GetNbrNId(i));
      if (max_out_deg < CurrI.GetOutDeg()) {
        max_out_deg = CurrI.GetOutDeg();
      }
      if (min_out_deg > CurrI.GetOutDeg()) {
        min_out_deg = CurrI.GetOutDeg();
      }
      sum_out_deg += int64(CurrI.GetOutDeg());
      MemNeeded += (int64(CurrI.GetOutDeg()))*(int64(sizeof(TInt) + sizeof(TFlt)));
      int64 tmp = (int64(CurrI.GetOutDeg()))*(int64(sizeof(TInt) + sizeof(TFlt)));
      // fprintf(fout, "%lld\t%lld\n", int64(CurrI.GetOutDeg()), tmp);
      if (tmp > 200000) {
        count_200k += 1;
        printf("tmp = %lld\n", tmp);
      }
    }
  }
  // fclose(fout);

  printf("min_out_deg= %d\n", min_out_deg);
  printf("max_out_deg= %d\n", max_out_deg);
  printf("sum_out_deg= %lld\n", sum_out_deg);
  printf("count_200k = %lld\n", count_200k);
  printf("count_node, count_edge = %d\t%d", count_node, count_edge);
  return MemNeeded;
}

//Simulates a random walk
void SimulateWalk(PWNet& InNet, int64 StartNId, int& WalkLen, TRnd& Rnd, TIntV& WalkV) {
  WalkV.Add(StartNId);
  if (WalkLen == 1) { return; }
  if (InNet->GetNI(StartNId).GetOutDeg() == 0) { return; }
  WalkV.Add(InNet->GetNI(StartNId).GetNbrNId(Rnd.GetUniDevInt(InNet->GetNI(StartNId).GetOutDeg())));
  while (WalkV.Len() < WalkLen) {
    int64 Dst = WalkV.Last();
    int64 Src = WalkV.LastLast();
    if (InNet->GetNI(Dst).GetOutDeg() == 0) { return; }
    int64 Next = AliasDrawInt(InNet->GetNDat(Dst).GetDat(Src),Rnd);
    WalkV.Add(InNet->GetNI(Dst).GetNbrNId(Next));
  }
}
