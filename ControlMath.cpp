
/****************************************************************************
 *
 *   Copyright (C) 2018 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file ControlMath.cpp
 */
 
/**
 * @pronouncement:ControlMath libriry is transplanted from the px4 firmaware to the MDK5 project. Obey to the license of BSD
 * @author snowpang 18516614106@163.com. All rights reserved.
 * @data 2018.9.26
 * @version ControlMath lib V1.0.0
 */

#include "ControlMath.hpp"
#include "mathlib.h"

namespace ControlMath
{
matrix::Quatf thrustToAttitude(const matrix::Vector3f &thr_sp, const float yaw_sp)
{
	
	float yaw_body_sp = yaw_sp;

	// desired body_z axis = -normalize(thrust_vector)
	matrix::Vector3f body_x, body_y, body_z;

	if (thr_sp.length() > 0.00001f) {
		body_z = -thr_sp.normalized();

	} else {
		// no thrust, set Z axis to safe value
		body_z.zero();
		body_z(2) = 1.0f;
	}

	// vector of desired yaw direction in XY plane, rotated by PI/2
	matrix::Vector3f y_C(-sinf(yaw_body_sp), cosf(yaw_body_sp), 0.0f);

	if (fabsf(body_z(2)) > 0.000001f) {
		// desired body_x axis, orthogonal to body_z
		body_x = y_C % body_z;

		// keep nose to front while inverted upside down
		if (body_z(2) < 0.0f) {
			body_x = -body_x;
		}

		body_x.normalize();

	} else {
		// desired thrust is in XY plane, set X downside to construct correct matrix,
		// but yaw component will not be used actually
		body_x.zero();
		body_x(2) = 1.0f;
	}

	// desired body_y axis
	body_y = body_z % body_x;

	matrix::Dcmf R_sp;

	// fill rotation matrix
	for (int i = 0; i < 3; i++) {
		R_sp(i, 0) = body_x(i);
		R_sp(i, 1) = body_y(i);
		R_sp(i, 2) = body_z(i);
	}

	//copy quaternion setpoint to attitude setpoint topic
	matrix::Quatf q_sp = R_sp;
//	q_sp.copyTo(att_sp.q_d);
//	att_sp.q_d_valid = true;

//	// calculate euler angles, for logging only, must not be used for control
//	matrix::Eulerf euler = R_sp;
//	att_sp.roll_body = euler(0);
//	att_sp.pitch_body = euler(1);
//	att_sp.thrust = thr_sp.length();

	return q_sp;
}

matrix::Vector2f constrainXY(const matrix::Vector2f &v0, const matrix::Vector2f &v1, const float &max)
{
	if (matrix::Vector2f(v0 + v1).norm() <= max) {
		// vector does not exceed maximum magnitude
		return v0 + v1;

	} else if (v0.length() >= max) {
		// the magnitude along v0, which has priority, already exceeds maximum.
		return v0.normalized() * max;

	} else if (fabsf(matrix::Vector2f(v1 - v0).norm()) < 0.001f) {
		// the two vectors are equal
		return v0.normalized() * max;

	} else if (v0.length() < 0.001f) {
		// the first vector is 0.
		return v1.normalized() * max;

	} else {
		// vf = final vector with ||vf|| <= max
		// s = scaling factor
		// u1 = unit of v1
		// vf = v0 + v1 = v0 + s * u1
		// constraint: ||vf|| <= max
		//
		// solve for s: ||vf|| = ||v0 + s * u1|| <= max
		//
		// Derivation:
		// For simplicity, replace v0 -> v, u1 -> u
		// 				   		   v0(0/1/2) -> v0/1/2
		// 				   		   u1(0/1/2) -> u0/1/2
		//
		// ||v + s * u||^2 = (v0+s*u0)^2+(v1+s*u1)^2+(v1+s*u1)^2 = max^2
		// v0^2+2*s*u0*v0+s^2*u0^2 + v1^2+2*s*u1*v1+s^2*u1^2 + v2^2+2*s*u2*v2+s^2*u2^2 = max^2
		// s^2*(u0^2+u1^2+u2^2) + s*2*(u0*v0+u1*v1+u2*v2) + (v0^2+v1^2+v2^2-max^2) = 0
		//
		// quadratic equation:
		// -> s^2*a + s*b + c = 0 with solution: s1/2 = (-b +- sqrt(b^2 - 4*a*c))/(2*a)
		//
		// b = 2 * u.dot(v)
		// a = 1 (because u is normalized)
		// c = (v0^2+v1^2+v2^2-max^2) = -max^2 + ||v||^2
		//
		// sqrt(b^2 - 4*a*c) =
		// 		sqrt(4*u.dot(v)^2 - 4*(||v||^2 - max^2)) = 2*sqrt(u.dot(v)^2 +- (||v||^2 -max^2))
		//
		// s1/2 = ( -2*u.dot(v) +- 2*sqrt(u.dot(v)^2 - (||v||^2 -max^2)) / 2
		//      =  -u.dot(v) +- sqrt(u.dot(v)^2 - (||v||^2 -max^2))
		// m = u.dot(v)
		// s = -m + sqrt(m^2 - c)
		//
		//
		//
		// notes:
		// 	- s (=scaling factor) needs to be positive
		// 	- (max - ||v||) always larger than zero, otherwise it never entered this if-statement

		matrix::Vector2f u1 = v1.normalized();
		float m = u1.dot(v0);
		float c = v0.length() * v0.length() - max * max;
		float s = -m + sqrtf(m * m - c);
		return v0 + u1 * s;
	}
}

matrix::Vector3f attitude_control(matrix::Quatf quat_desire, matrix::Quatf quat_current, float yaw_w)
{
	/* 1. get estimated and desired vehicle attitude */
	matrix::Quatf q = quat_current;
	matrix::Quatf qd = quat_desire;
	
	/* ensure input quaternions are exactly normalized because acosf(1.00001) == NaN */
	q.normalize();
	qd.normalize();
	
	/* 3.calculate reduced desired attitude neglecting vehicle's yaw to prioritize roll and pitch */
	matrix::Vector3f e_z = q.dcm_z();
	matrix::Vector3f e_z_d = qd.dcm_z();
	matrix::Quatf qd_red(e_z, e_z_d);
	
	if (abs(qd_red(1)) > (1.f - 1e-5f) || abs(qd_red(2)) > (1.f - 1e-5f)) { 
		/* In the infinitesimal corner case where the vehicle and thrust have the completely opposite direction,
		 * full attitude control anyways generates no yaw input and directly takes the combination of
		 * roll and pitch leading to the correct desired yaw. Ignoring this case would still be totally safe and stable. */
		qd_red = qd;

	} else {
		/* transform rotation from current to desired thrust vector into a world frame reduced desired attitude */
		qd_red *= q;
	}
	
	/* 4.mix full and reduced desired attitude */
	matrix::Quatf q_mix = qd_red.inversed() * qd;
	/* find the shortest rotate route (+-pi) */
	q_mix *= math::signNoZero(q_mix(0)); 
	/* catch numerical problems with the domain of acosf and asinf */
	q_mix(0) = math::constrain(q_mix(0), -1.f, 1.f);
	q_mix(3) = math::constrain(q_mix(3), -1.f, 1.f);
	/* 5. recalculate the yaw control output in reduced by yaw weight*/
	qd = qd_red * matrix::Quatf(cosf(yaw_w * acosf(q_mix(0))), 0, 0, sinf(yaw_w * asinf(q_mix(3))));

	/*6. quaternion attitude control law, qe is rotation from q to qd */
	matrix::Quatf qe = q.inversed() * qd;

	/*7. using sin(alpha/2) scaled rotation axis as attitude error (see quaternion definition by axis angle)
	 * also taking care of the antipodal unit quaternion ambiguity */
	matrix::Vector3f eq = 2.f * math::signNoZero(qe(0)) * qe.imag();
	return eq;
}

}
