#ifndef _Record_H
#define _Record_H

#include<Vehicle_State.h>
#include<Controller_State.h>
#include<Guidance_State.h>

void Record(const VehicleState& Vehicle_State, const ControllerState& Controller_State, const GuidanceState Guidance_State, FILE *File_Vehicle, FILE *File_Control);

#endif