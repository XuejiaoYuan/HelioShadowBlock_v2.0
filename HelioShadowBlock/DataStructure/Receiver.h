#pragma once
//
// Created by Amber on 2018/4/3.
//
//  Receiver
//  Define the receiver in the solar system.
//

#ifndef HELIOSHADOW_RECEIVER_H
#define HELIOSHADOW_RECEIVER_H

//#include "../Common/global_function.h"
#include "../Common/CommonFunc.h"

typedef enum {
	RectangularRecvType, CylinderRecvType, CircularTruncatedConeRecvType
}ReceiverType;

class Receiver {
public:
	Receiver(const ReceiverType&_recv_type) {
		recv_type = _recv_type;
		recv_pos = Vector3f(0, 0, 0);
		recv_size = Vector3f(0, 0, 0);
		recv_normal = Vector3f(0, 0, 0);
		recv_face = 0;
	}

	ReceiverType recv_type;             //Receiver's type
	Vector3f recv_pos;                      //The position of receiver's center
	Vector3f recv_size;                     //Receiver's size
	Vector3f recv_normal;                   //The normal of the receiver's face
	int recv_face;                      //The index of the receiver's face

};

class RectangularRecv :public Receiver {
public:
	RectangularRecv() :Receiver(RectangularRecvType) {};
};

class CylinderRecv :public Receiver {
public:
	CylinderRecv() :Receiver(CylinderRecvType) {};
};

class CircularTruncatedConeRecv :public Receiver {
public:
	CircularTruncatedConeRecv() :Receiver(CircularTruncatedConeRecvType) {};
};

class ReceiverCreator {
public:
	Receiver* getReceiver(const ReceiverType& recv_type) {
		switch (recv_type) {
		case RectangularRecvType:
			return new RectangularRecv();
		case CylinderRecvType:
			return new CylinderRecv();
		case CircularTruncatedConeRecvType:
			return new CircularTruncatedConeRecv();
		default:
			return new RectangularRecv();
		}
	}
};
#endif //HELIOSHADOW_RECEIVER_H
