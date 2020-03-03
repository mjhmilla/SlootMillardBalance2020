return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {25.0},
	 subject_height = {1.65},
	 subject_weight = {60.00},
	 subject_gender = {"female"},
	 subject_pelvisWidth = {0.2100},
	 subject_hipCenterWidth = {0.1600},
	 subject_shoulderCenterWidth = {0.2460},
	 subject_heelAnkleXOffset = {0.0990},
	 subject_heelAnkleZOffset = {0.0950},
	 subject_shoulderNeckZOffset = {0.0590},
	 subject_footWidth = {0.1200},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Pelvis_L", body = "Pelvis", point = {0.000000, 0.172608, 0.120826,},},
	{name = "Pelvis_R", body = "Pelvis", point = {0.000000, -0.172608, 0.120826,},},
	{name = "Pelvis_Back", body = "Pelvis", point = {-0.129456, 0.000000, 0.120826,},},
	{name = "Pelvis_Front", body = "Pelvis", point = {0.129456, 0.000000, 0.120826,},},
	{name = "Thigh_R", body = "Thigh_R", point = {0.070089, -0.000000, -0.262835,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.099000, 0.060000, -0.095000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.099000, -0.060000, -0.095000,},},
	{name = "Toe_R", body = "Foot_R", point = {0.063836, 0.000000, -0.095000,},},
	{name = "Thigh_L", body = "Thigh_L", point = {0.070089, 0.000000, -0.262835,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.099000, -0.060000, -0.095000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.099000, 0.060000, -0.095000,},},
	{name = "Toe_L", body = "Foot_L", point = {0.063836, 0.000000, -0.095000,},},
	{name = "UpperTrunk_Front", body = "UpperTrunk", point = {0.086732, 0.000000, 0.065049,},},
	{name = "UpperTrunk_Back", body = "UpperTrunk", point = {-0.097573, 0.000000, 0.086732,},},
	{name = "ProximalMetacarpal_Medial_R", body = "Hand_R", point = {-0.032353, 0.024265, -0.032353,},},
	{name = "ProximalMetacarpal_Lateral_R", body = "Hand_R", point = {0.032353, 0.024265, -0.032353,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.032353, 0.024265, -0.097060,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.032353, 0.024265, -0.097060,},},
	{name = "ProximalMetacarpal_Medial_L", body = "Hand_L", point = {-0.032353, -0.024265, -0.032353,},},
	{name = "ProximalMetacarpal_Lateral_L", body = "Hand_L", point = {0.032353, -0.024265, -0.032353,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.032353, -0.024265, -0.097060,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.032353, -0.024265, -0.097060,},},
},
constraint_sets = {
},
frames = {
	{name   = "Pelvis",
	parent = "ROOT",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 7.482000,
		com = 	{ 0.000000, 0.000000, 0.087685,},
		inertia = 
			{{ 0.041794, 0.000000, 0.000000,},
			{ 0.000000, 0.036024, 0.000000,},
			{ 0.000000, 0.000000, 0.043945,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.129456, 0.086304, 0.103565,},
		R_IAS = 	{ 0.129456, -0.086304, 0.103565,},
		L_IPS = 	{ -0.120826, 0.120826, 0.120826,},
		R_IPS = 	{ -0.120826, -0.120826, 0.120826,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.210000, 0.262500, 0.258912,},
		mesh_center = 	{ 0.000000, 0.000000, 0.043152,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.080000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 8.868000,
		com = 	{ 0.000000, 0.000000, -0.126581,},
		inertia = 
			{{ 0.148294, 0.000000, 0.000000,},
			{ 0.000000, 0.144302, 0.000000,},
			{ 0.000000, 0.000000, 0.028582,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.070089, -0.115647,},
		R_FLE = 	{ 0.000000, -0.052567, -0.350447,},
		R_FME = 	{ 0.000000, 0.052567, -0.350447,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.175223, 0.140179, 0.420536,},
		mesh_center = 	{ 0.000000, 0.000000, -0.175223,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.350447,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.886000,
		com = 	{ 0.000000, 0.000000, -0.181527,},
		inertia = 
			{{ 0.035795, 0.000000, 0.000000,},
			{ 0.000000, 0.034731, 0.000000,},
			{ 0.000000, 0.000000, 0.004250,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.083422, -0.020856,},
		R_TTC = 	{ 0.020856, 0.000000, -0.062567,},
		R_FAL = 	{ 0.000000, -0.062567, -0.417112,},
		R_TAM = 	{ 0.000000, 0.062567, -0.417112,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.125134, 0.125134, 0.500535,},
		mesh_center = 	{ 0.000000, 0.000000, -0.208556,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.417112,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.774000,
		com = 	{ -0.011850, 0.000000, -0.047500,},
		inertia = 
			{{ 0.003262, 0.000000, 0.000000,},
			{ 0.000000, 0.002840, 0.000000,},
			{ 0.000000, 0.000000, 0.000705,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.151981, 0.021712, -0.032567,},
		R_FM2 = 	{ 0.151981, 0.000000, -0.032567,},
		R_FM5 = 	{ 0.151981, -0.021712, -0.032567,},
		R_FCC = 	{ -0.065135, 0.000000, -0.032567,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.217115, 0.120000, 0.095000,},
		mesh_center = 	{ 0.009558, 0.000000, -0.047500,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.080000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 8.868000,
		com = 	{ 0.000000, 0.000000, -0.126581,},
		inertia = 
			{{ 0.148294, 0.000000, 0.000000,},
			{ 0.000000, 0.144302, 0.000000,},
			{ 0.000000, 0.000000, 0.028582,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.070089, -0.115647,},
		L_FLE = 	{ 0.000000, 0.052567, -0.350447,},
		L_FME = 	{ 0.000000, -0.052567, -0.350447,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.175223, 0.140179, 0.420536,},
		mesh_center = 	{ 0.000000, 0.000000, -0.175223,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.350447,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.886000,
		com = 	{ 0.000000, 0.000000, -0.181527,},
		inertia = 
			{{ 0.035795, 0.000000, 0.000000,},
			{ 0.000000, 0.034731, 0.000000,},
			{ 0.000000, 0.000000, 0.004250,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.083422, -0.020856,},
		L_TTC = 	{ 0.020856, 0.000000, -0.062567,},
		L_FAL = 	{ 0.000000, 0.062567, -0.417112,},
		L_TAM = 	{ 0.000000, -0.062567, -0.417112,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.125134, 0.125134, 0.500535,},
		mesh_center = 	{ 0.000000, 0.000000, -0.208556,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.417112,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.774000,
		com = 	{ -0.011850, 0.000000, -0.047500,},
		inertia = 
			{{ 0.003262, 0.000000, 0.000000,},
			{ 0.000000, 0.002840, 0.000000,},
			{ 0.000000, 0.000000, 0.000705,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.151981, -0.021712, -0.032567,},
		L_FM2 = 	{ 0.151981, 0.000000, -0.032567,},
		L_FM5 = 	{ 0.151981, 0.021712, -0.032567,},
		L_FCC = 	{ -0.065135, 0.000000, -0.032567,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.217115, 0.120000, 0.095000,},
		mesh_center = 	{ 0.009558, 0.000000, -0.047500,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.172608,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 8.790000,
		com = 	{ 0.000000, 0.000000, 0.107149,},
		inertia = 
			{{ 0.062822, 0.000000, 0.000000,},
			{ 0.000000, 0.041990, 0.000000,},
			{ 0.000000, 0.000000, 0.057707,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.097621, 0.000000, 0.097621,},
		MAI = 	{ -0.097621, 0.000000, 0.097621,},
		LV1 = 	{ -0.097621, 0.000000, 0.000000,},
		LV3 = 	{ -0.097621, 0.000000, -0.048811,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.189000, 0.231000, 0.292863,},
		mesh_center = 	{ 0.000000, 0.000000, 0.048811,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.195242,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 9.270000,
		com = 	{ 0.000000, 0.000000, 0.107331,},
		inertia = 
			{{ 0.094643, 0.000000, 0.000000,},
			{ 0.000000, 0.042971, 0.000000,},
			{ 0.000000, 0.000000, 0.087864,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.086732, 0.000000, 0.173464,},
		SJN = 	{ 0.091069, 0.000000, 0.130098,},
		TV2 = 	{ -0.086732, 0.000000, 0.130098,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.161438, 0.307500, 0.227671,},
		mesh_center = 	{ 0.000000, 0.000000, 0.086732,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.216830,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 4.008000,
		com = 	{ 0.000000, 0.000000, 0.119565,},
		inertia = 
			{{ 0.015811, 0.000000, 0.000000,},
			{ 0.000000, 0.018735, 0.000000,},
			{ 0.000000, 0.000000, 0.014665,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.069528, 0.139056,},
		R_HEAD = 	{ 0.000000, -0.069528, 0.139056,},
		SGL    = 	{ 0.069528, 0.000000, 0.139056,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.185409, 0.185409, 0.243349,},
		mesh_center = 	{ 0.000000, 0.000000, 0.092704,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.123000, 0.157830,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.530000,
		com = 	{ 0.000000, 0.000000, -0.150538,},
		inertia = 
			{{ 0.008093, 0.000000, 0.000000,},
			{ 0.000000, 0.007079, 0.000000,},
			{ 0.000000, 0.000000, 0.002294,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.026162,},
		R_HUM = 	{ 0.000000, -0.039243, -0.175287,},
		R_HLE = 	{ 0.000000, -0.039243, -0.261622,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.130811, 0.104649, 0.287785,},
		mesh_center = 	{ 0.000000, 0.000000, -0.130811,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.261622,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.828000,
		com = 	{ 0.000000, 0.000000, -0.114591,},
		inertia = 
			{{ 0.003563, 0.000000, 0.000000,},
			{ 0.000000, 0.003455, 0.000000,},
			{ 0.000000, 0.000000, 0.000462,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ 0.050270, 0.000000, -0.251352,},
		R_RSP = 	{ -0.050270, 0.000000, -0.251352,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.100541, 0.075405, 0.276487,},
		mesh_center = 	{ 0.000000, 0.000000, -0.125676,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.251352,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.336000,
		com = 	{ 0.000000, 0.000000, -0.055437,},
		inertia = 
			{{ 0.000523, 0.000000, 0.000000,},
			{ 0.000000, 0.000380, 0.000000,},
			{ 0.000000, 0.000000, 0.000209,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.129413,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.113237, 0.048530, 0.161767,},
		mesh_center = 	{ 0.000000, 0.000000, -0.080883,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.123000, 0.157830,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.530000,
		com = 	{ 0.000000, 0.000000, -0.150538,},
		inertia = 
			{{ 0.008093, 0.000000, 0.000000,},
			{ 0.000000, 0.007079, 0.000000,},
			{ 0.000000, 0.000000, 0.002294,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.026162,},
		L_HUM = 	{ 0.000000, 0.039243, -0.130811,},
		L_HLE = 	{ 0.000000, 0.039243, -0.261622,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.130811, 0.104649, 0.287785,},
		mesh_center = 	{ 0.000000, -0.000000, -0.130811,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.261622,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.828000,
		com = 	{ 0.000000, 0.000000, -0.114591,},
		inertia = 
			{{ 0.003563, 0.000000, 0.000000,},
			{ 0.000000, 0.003455, 0.000000,},
			{ 0.000000, 0.000000, 0.000462,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.050270, 0.000000, -0.251352,},
		L_RSP = 	{ 0.050270, 0.000000, -0.251352,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.100541, 0.075405, 0.276487,},
		mesh_center = 	{ 0.000000, 0.000000, -0.125676,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.251352,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.336000,
		com = 	{ 0.000000, 0.000000, -0.055437,},
		inertia = 
			{{ 0.000523, 0.000000, 0.000000,},
			{ 0.000000, 0.000380, 0.000000,},
			{ 0.000000, 0.000000, 0.000209,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.129413,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.113237, 0.048530, 0.161767,},
		mesh_center = 	{ 0.000000, 0.000000, -0.080883,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}