#ifndef _Vehicle_State_Update_H
#define _Vehicle_State_Update_H

#include<Global.h>

struct VehicleState Vehicle_State_Update(struct VehicleState Vehicle_State,struct ControllerState Controller_State,struct VehiclePara Vehicle_Para,double h);

#endif