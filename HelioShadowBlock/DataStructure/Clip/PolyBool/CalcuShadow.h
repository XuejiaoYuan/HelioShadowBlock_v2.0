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


struct DLLEXPORT Polyg   //ƽ������
{
	Point2D p[15];                
	int vertNum;
}; 

struct DLLEXPORT EffectiveRegion
{
	int PolNum;        //����������θ���
    Polyg pol[5];      //�໥�����ƽ������
	float LitArea;     //���������
};

        
DLLEXPORT EffectiveRegion ShadowBlockCalcuCPU( const float a_SunDir[ 3 ], const float aFocus_pos[3], int ReflNum, int CurMirIndex,  CHeliostat *a_pReflector); //���ع���������μ����������

DLLEXPORT bool DetectShadowMirror( CHeliostat testRefl, CHeliostat curRefl, float hitPos[4][3], float RayDir[3] );       //true������ڵ�

//DLLEXPORT float SpillageLoss( const float recv_pos[3], const float recv_Goem[3], const float (*Recv_rotate_Matrix)[4], const CReflector cur_Reflector , EffectiveRegion LightShape);        //�������������Ľ���Ч�ʣ��������+����ЧӦ��

DLLEXPORT float SpillageLoss( CReceiver* mReceiver, const CHeliostat cur_Reflector, EffectiveRegion LightShape);