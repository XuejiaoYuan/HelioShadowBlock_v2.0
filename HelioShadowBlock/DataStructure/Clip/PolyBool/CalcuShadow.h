#pragma once
#ifndef DLLEXPORT
#define DLLEXPORT _declspec(dllexport)
#endif

#include "../../DataCollection/SolarScene.h"
#include "../../DataCollection/MirrorPlane.h"
#include "../culling/Struct.h"


#include <math.h>
#include <cstring>
using namespace::std;


typedef struct DLLEXPORT Point2D{double x,y;}Vector2D,Point2D;
#define N_EDGE 4


struct DLLEXPORT Polyg   //平面多边形
{
	Point2D p[15];                
	int vertNum;
}; 

struct DLLEXPORT EffectiveRegion
{
	int PolNum;        //光亮区多边形个数
    Polyg pol[5];      //相互分离的平面多边形
	float LitArea;     //光亮区面积
};

        
DLLEXPORT EffectiveRegion ShadowBlockCalcuCPU( const float a_SunDir[ 3 ], const float aFocus_pos[3], int ReflNum, int CurMirIndex,  CHeliostat *a_pReflector); //返回光亮区多边形及光亮区面积

DLLEXPORT bool DetectShadowMirror( CHeliostat testRefl, CHeliostat curRefl, float hitPos[4][3], float RayDir[3] );       //true代表会遮挡

//DLLEXPORT float SpillageLoss( const float recv_pos[3], const float recv_Goem[3], const float (*Recv_rotate_Matrix)[4], const CReflector cur_Reflector , EffectiveRegion LightShape);        //计算接收器表面的接收效率（能量溢出+余弦效应）

DLLEXPORT float SpillageLoss( CReceiver* mReceiver, const CHeliostat cur_Reflector, EffectiveRegion LightShape);