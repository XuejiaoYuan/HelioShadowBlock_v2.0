#include "CalcuShadow.h"

#include <time.h>  //��ʱ

#include "PolyBool/poly.h"
#include "PolyBool/poly_io.h"
#include "PolyBool/nclip.h"

#define PI 3.1415926535


// inline void clear(PolyPList &l)
// {
// 	PolyPListIter	i(l);
// 	while(i())
// 		delete i.val();
// }

inline void CoordTrans( float orig[3],  float xaxis[3],  float yaxis[3],  float zaxis[3], float (*in_hitPos)[3], Point out_point2D[4])
{
	for(int vertCount =0; vertCount<4; vertCount++)
	{
		for(int dim =0; dim<3; dim++ )
			in_hitPos[vertCount][dim] -=orig[dim];          //ƽ��
		//��ת
		float xcoord = DOT(in_hitPos[vertCount], xaxis);
		float ycoord = DOT(in_hitPos[vertCount], yaxis);
		out_point2D[vertCount] = Point(xcoord, ycoord);
	}
}

DLLEXPORT bool DetectShadowMirror( CHeliostat testRefl, CHeliostat curRefl, float hitPos[4][3], float RayDir[3] )       //true����testRefl���ڵ�curRefl
{
	int i;
	const float *pos, *geo;
	float (*mRotate)[4];
	float rayOrig[4][3];       //���쾵���ĸ����㷢���Ĺ���

	pos = curRefl.GetPos();
	geo = curRefl.GetGeo();
	mRotate = curRefl.GetRotateMatrix();

	MirrorPlane curMirror = MirrorPlane( pos, geo[0], geo[1], mRotate);

	pos = testRefl.GetPos();
	geo = testRefl.GetGeo();
	mRotate = testRefl.GetRotateMatrix();

	MirrorPlane testMirror = MirrorPlane(pos, geo[0], geo[1], mRotate);

	for( i = 0; i < 3; i++)
	{
		rayOrig[0][i] = testMirror.getPos()[i] - 0.5 * testMirror.getW() * testMirror.getXaxis()[i] +  0.5 * testMirror.getH() * testMirror.getYaxis()[i];
		rayOrig[1][i] = testMirror.getPos()[i] + 0.5 * testMirror.getW() * testMirror.getXaxis()[i] +  0.5 * testMirror.getH() * testMirror.getYaxis()[i];
		rayOrig[2][i] = testMirror.getPos()[i] + 0.5 * testMirror.getW() * testMirror.getXaxis()[i] -  0.5 * testMirror.getH() * testMirror.getYaxis()[i];
		rayOrig[3][i] = testMirror.getPos()[i] - 0.5 * testMirror.getW() * testMirror.getXaxis()[i] -  0.5 * testMirror.getH() * testMirror.getYaxis()[i];
	}

	bool IsShadow = false;              //      ��ǲ��Ծ����Ƿ���ڵ���ǰ����

	//ray-mirror hit test  
	for(i=0; i<4; i++)
		if( curMirror.RayHit(rayOrig[i], RayDir, hitPos[i] ) )         
		{
			if( !IsShadow )
				IsShadow = true;
		}
      
	return IsShadow;
}


EffectiveRegion ShadowBlockCalcuCPU( const float a_SunDir[ 3 ], const float aFocus_pos[3], int ReflNum,  int CurMirIndex, CHeliostat *a_pReflector) //���ع����������      //20130826���
{
	float invSun[3] = { -a_SunDir[0], -a_SunDir[1], -a_SunDir[2] };            //̫��������

	//  debug
// //	printf("the sun direction in ShadowCalc model: %f,%f,%f", vSun.x, vSun.y, vSun.z);   
// 
// // 	int i;           
// // 
//      
// 
	const float *pos, *geo;
	float (*mRotate)[4];
	float rayOrig[4][3];       //���쾵���ĸ����㷢���Ĺ���
 	float hitPos[4][3];         //�����ĸ������ͶӰ����3D     

	MirrorPlane curMirror;

	EffectiveRegion CurMirrorRes;
	CurMirrorRes.LitArea = 0.0;    //��ʼ��   
  
	PolyPList LightPart;          //�洢�����������
    
	CHeliostat curRefl = a_pReflector[CurMirIndex];  
               //creat the virtual sun frustum 
	
	pos = curRefl.GetPos();
	geo = curRefl.GetGeo();
	mRotate = curRefl.GetRotateMatrix();

	curMirror = MirrorPlane( pos, geo[0], geo[1], mRotate);

	Vector3 reflDir = Vector3(aFocus_pos[0] - pos[0],aFocus_pos[1] - pos[1], aFocus_pos[2] - pos[2] );
	reflDir = reflDir.NormV3();          //��һ��
	float inv_ref[3] = { -reflDir.x, -reflDir.y, -reflDir.z };
       
	//�����ʼ��������
	Point p[4];
	p[0] = Point(geo[0]/2.0, geo[1]/2.0), p[1] = Point(geo[0]/2.0, -geo[1]/2.0),   
		p[2] = Point(-geo[0]/2.0, -geo[1]/2.0), p[3] = Point(-geo[0]/2.0, geo[1]/2.0);
	Poly origPoly = Poly( p[0] );
	for(int vert =1; vert < 4; vert++)
		origPoly.add(p[vert]);
	LightPart.add(&origPoly);  

// 	vector<int> SBIndex = a_pReflector[CurMirIndex].GetSBVect();
// 	cout<<"the shadow block list size: "<<SBIndex.size()<<endl;

	vector<int> SBIndex;               //  ɾ������ӰԤ��������û������Ӱ�ڵ������Ժ��õ��������         20140618

//	for(vector<int>::iterator it = SBIndex.begin(); it!=SBIndex.end(); it++)
	for (int MirIndex = 0; MirIndex < ReflNum; MirIndex++)	
	{
//		int MirIndex = *it;
		if( MirIndex == CurMirIndex )     
			continue;

		if( DetectShadowMirror(a_pReflector[MirIndex], a_pReflector[CurMirIndex], hitPos, invSun))       //���Ծ��治�������Ӱ  shadow����
		{
//			printf("the shadow mirror index %d\n", MirIndex);
			//����仯
			Point ShadowCorner[4];
			CoordTrans(curMirror.getPos(), curMirror.getXaxis(), curMirror.getYaxis(),curMirror.getNormal(), hitPos, ShadowCorner);  //����任

			//�ü�
			Poly ShadowShape = Poly(ShadowCorner[0]);
			for(int vert1 =1; vert1 < 4; vert1++ )
				ShadowShape.add(ShadowCorner[vert1]);
			PolyPList e_min_d,d_min_e,d_and_e;
			PolyPListIter iter(LightPart);
			while(iter())
				clip_poly(*iter.val(), ShadowShape, e_min_d, d_min_e, d_and_e);    //kernel
			// 		clear(d_min_e);
			// 		clear(d_and_e); 
			//���ԭ��������״����
			PolyPListIter iter1(LightPart);
			while(iter1())
				LightPart.del(iter1.val());
			PolyPListIter iter2(LightPart);
			while(iter2())
				LightPart.del(iter2.val());

			//��������״������Ӳü�������
			PolyPListIter iter3(e_min_d);
			while(iter3())
				LightPart.add(iter3.val());
		}

		if (DetectShadowMirror(a_pReflector[MirIndex], a_pReflector[CurMirIndex], hitPos, inv_ref))
		{
//			printf("the block mirror index %d\n", MirIndex);
			//����仯
			Point BlockCorner[4];
			CoordTrans(curMirror.getPos(), curMirror.getXaxis(), curMirror.getYaxis(),curMirror.getNormal(), hitPos, BlockCorner);  //����任

			//�ü�
			Poly BlockShape = Poly(BlockCorner[0]);
			for(int vert1 =1; vert1 < 4; vert1++ )
				BlockShape.add(BlockCorner[vert1]);

			PolyPList e_min_d,d_min_e,d_and_e;
			PolyPListIter iter4(LightPart);
			while(iter4())
				clip_poly(*iter4.val(), BlockShape, e_min_d, d_min_e, d_and_e);    //kernel
			// 		clear(d_min_e);
			// 		clear(d_and_e); 
			//���ԭ��������״����
			PolyPListIter iter5(LightPart);
			while(iter5())
				LightPart.del(iter5.val());
			PolyPListIter iter6(LightPart);
			while(iter6())
				LightPart.del(iter6.val());

			//��������״������Ӳü�������
			PolyPListIter iter7(e_min_d);
			while(iter7())
				LightPart.add(iter7.val());
		}
		
	}
	

	PolyPListIter iter(LightPart);
	int poly_count = 0;
	while (iter())
	{
		Poly temp = Poly(*iter.val());
		ConstPolyIter iter(temp);
		int vert = 0;
		while(iter())
		{
			CurMirrorRes.pol[poly_count].p[vert].x = iter.point().x();
			CurMirrorRes.pol[poly_count].p[vert].y = iter.point().y();
			vert++;
		}
		CurMirrorRes.pol[poly_count].vertNum = vert;
		poly_count++;
		CurMirrorRes.LitArea += temp.area();        //���ú����������
	}
	CurMirrorRes.PolNum = poly_count;

	return CurMirrorRes;
}   



float SpillageLoss(  CReceiver* mReceiver, const CHeliostat cur_Reflector, EffectiveRegion LightShape)     //������ЧͶӰ�ı���
{

// 	float ratio = 0;
// 	float EffectArea = 0.0, ProjectArea = 0.0;       //��Ч���ཻ���, ����ֻ������Ч�������������״
// 
// 	//��õ�ǰ���վ�����
// 	CHeliostat mReflector = cur_Reflector;
// 	const float* Refl_pos = mReflector.GetPos();
// 	float (*Refl_rotate_Matrix)[4] = mReflector.GetRotateMatrix();
// 
// 	Vector3 Refl_pos_v3 = Vector3(Refl_pos[0], Refl_pos[1], Refl_pos[2]);
// 	Vector3 Refl_Xaxis_v3 = Vector3(Refl_rotate_Matrix[0][0],Refl_rotate_Matrix[0][1], Refl_rotate_Matrix[0][2] );
// 	Vector3 Refl_Yaxis_v3 = Vector3(Refl_rotate_Matrix[2][0],Refl_rotate_Matrix[2][1], Refl_rotate_Matrix[2][2] );
// 	Vector3 Refl_Norm_v3 = Vector3(Refl_rotate_Matrix[1][0],Refl_rotate_Matrix[1][1], Refl_rotate_Matrix[1][2]);    
// 
// 	//���㷴�䷽��
// 
// 	Vector3 Focus_point = mReceiver->GetFocusPoint();
// 	Vector3 ReflDir = Focus_point - Refl_pos_v3;
// 	ReflDir = ReflDir.NormV3();
// 
// 	ReceivePanel* mPanel = mReceiver->GetPanel();
// 
// 	//�Խ������������ѭ��
//     for (int PanelIdex = 0; PanelIdex < mReceiver->GetPanelCount(); PanelIdex++)
//     {
// // 		if( !mPanel[PanelIdex].IsRayForward(ReflDir) )    //��ǰ�����ղ�������
// // 			continue;
// 		float pshapeW,pshapeH;
// 		mPanel[PanelIdex].GetPanelShape(pshapeW, pshapeH);   //��ά
// 
// 		Vector3* pCorner = mPanel[PanelIdex].GetPanelCorner();
// 
// 		Vector3 pNorm = mPanel[PanelIdex].GetPanelNorm();
// 
// 		//���챻�ü����������
// 		Point p[4];
// 		p[0] = Point(0.0, 0.0), p[1] = Point(pshapeW, 0),   
// 			p[2] = Point(pshapeW, pshapeH), p[3] = Point(0, pshapeH);
// 		Poly tailorPoly = Poly( p[0] );          
// 		for(int vert =1; vert < 4; vert++)
// 			tailorPoly.add(p[vert]);
// 
// 		for(int npol = 0; npol < LightShape.PolNum; npol++)    
// 		{
// 			Polyg m_poly = LightShape.pol[npol];
// 			Vector3* shapeCorner = new Vector3[m_poly.vertNum];
// 			Vector3* hitPos = new Vector3[m_poly.vertNum];          //�����������ط������ǰ���������Ľ���ռ�����
// 			Point* hitPoint = new Point[m_poly.vertNum];
// 			for(int vet = 0; vet < m_poly.vertNum; vet++)
// 			{
// 				shapeCorner[vet] = Refl_pos_v3 + m_poly.p[vet].x * Refl_Xaxis_v3 + m_poly.p[vet].y * Refl_Yaxis_v3; 
// 				
// 				//float k = ((recvPlanePoint * recvPlaneNorm) - (shapeCorner[vet]  * recvPlaneNorm)) * 1.0 /(reflDir_v3 * recvPlaneNorm);    //���㷴����߸�����ƽ��Ľ���  k = (p*n - o*n)/(d*n)
// 				float k = ((pCorner[0] * pNorm) - (shapeCorner[vet]  * pNorm)) * 1.0 /(ReflDir * pNorm);    //���㷴����߸�����ƽ��Ľ���  k = (p*n - o*n)/(d*n)
// 				hitPos[vet] = shapeCorner[vet] + k * ReflDir;
// 				hitPoint[vet] = Point( (hitPos[vet] - pCorner[0])* mPanel[PanelIdex].GetXaxis(), (hitPos[vet] - pCorner[0]) * mPanel[PanelIdex].GetZaxis());
// /*				int mmm = 0;*/
// 			}
// 
// 			Poly curPoly = Poly(hitPoint[0]);
// 			for(int vet = 1; vet < m_poly.vertNum; vet++)
// 				curPoly.add(hitPoint[vet]);
// 			ProjectArea += curPoly.area();
// 
// 			PolyPList e_min_d,d_min_e,d_and_e;
// 			clip_poly(tailorPoly, curPoly, e_min_d, d_min_e, d_and_e);    //kernel   �ü�  �����ཻ��������� ����Ч�Ľ������
// 
// 			//������Ӳü�������
// 			PolyPListIter iter(d_and_e);           //��ʱ�����ཻ���
// 			while(iter())
// 				EffectArea +=Poly(*iter.val()).area();
// 
// 			delete []shapeCorner;
// 			delete []hitPoint;
// 			delete []hitPos;
// 		}
// 
// 		ratio += (float)EffectArea / ProjectArea;
// 
//     }
// 
// 	return ratio;

}