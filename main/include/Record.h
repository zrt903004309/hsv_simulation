#ifndef _Record_H
#define _Record_H

#include<Vehicle_State.h>
#include<Controller_State.h>

void Record(VehicleState Vehicle_State, ControllerState Controller_State, FILE *File_Vehicle, FILE *File_Control);

#endif