return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {25.0},
	 subject_height = {1.68},
	 subject_weight = {63.00},
	 subject_gender = {"female"},
	 subject_pelvisWidth = {0.2340},
	 subject_hipCenterWidth = {0.1780},
	 subject_shoulderCenterWidth = {0.2010},
	 subject_heelAnkleXOffset = {0.0670},
	 subject_heelAnkleZOffset = {0.1000},
	 subject_shoulderNeckZOffset = {0.0310},
	 subject_footWidth = {0.1170},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Thigh_Proximal_R", body = "Thigh_R", point = {-0.071364, 0.000000, 0.035682,},},
	{name = "Thigh_Distal_R", body = "Thigh_R", point = {-0.071364, 0.000000, -0.089205,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.067000, 0.040950, -0.100000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.067000, -0.040950, -0.100000,},},
	{name = "ForeFoot_Medial_R", body = "Foot_R", point = {0.154063, 0.052650, -0.100000,},},
	{name = "ForeFoot_Lateral_R", body = "Foot_R", point = {0.154063, -0.052650, -0.100000,},},
	{name = "Thigh_Proximal_L", body = "Thigh_L", point = {-0.071364, 0.000000, 0.035682,},},
	{name = "Thigh_Distal_L", body = "Thigh_L", point = {-0.071364, 0.000000, -0.089205,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.067000, -0.040950, -0.100000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.067000, 0.040950, -0.100000,},},
	{name = "ForeFoot_Medial_L", body = "Foot_L", point = {0.154063, -0.052650, -0.100000,},},
	{name = "ForeFoot_Lateral_L", body = "Foot_L", point = {0.154063, 0.052650, -0.100000,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.032942, 0.008235, -0.082354,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.032942, 0.008235, -0.082354,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.032942, -0.008235, -0.082354,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.032942, -0.008235, -0.082354,},},
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
		mass   = 7.856100,
		com = 	{ 0.000000, 0.000000, 0.089279,},
		inertia = 
			{{ 0.045494, 0.000000, 0.000000,},
			{ 0.000000, 0.039213, 0.000000,},
			{ 0.000000, 0.000000, 0.047835,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.131810, 0.087873, 0.105448,},
		R_IAS = 	{ 0.131810, -0.087873, 0.105448,},
		L_IPS = 	{ -0.123022, 0.052724, 0.123022,},
		R_IPS = 	{ -0.123022, -0.052724, 0.123022,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.210600, 0.234000, 0.263620,},
		mesh_center = 	{ 0.000000, 0.000000, 0.043937,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.089000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 9.311400,
		com = 	{ 0.000000, 0.000000, -0.128883,},
		inertia = 
			{{ 0.161422, 0.000000, 0.000000,},
			{ 0.000000, 0.157077, 0.000000,},
			{ 0.000000, 0.000000, 0.031113,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.071364, -0.117750,},
		R_FLE = 	{ 0.000000, -0.053523, -0.356818,},
		R_FME = 	{ 0.000000, 0.053523, -0.356818,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.142727, 0.142727, 0.428182,},
		mesh_center = 	{ 0.000000, 0.000000, -0.178409,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.356818,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 3.030300,
		com = 	{ 0.000000, 0.000000, -0.184828,},
		inertia = 
			{{ 0.038964, 0.000000, 0.000000,},
			{ 0.000000, 0.037805, 0.000000,},
			{ 0.000000, 0.000000, 0.004626,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.084939, -0.021235,},
		R_TTC = 	{ 0.021235, 0.000000, -0.063704,},
		R_FAL = 	{ 0.000000, -0.063704, -0.424696,},
		R_TAM = 	{ 0.000000, 0.063704, -0.424696,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.127409, 0.127409, 0.509636,},
		mesh_center = 	{ 0.000000, 0.000000, -0.212348,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.424696,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.812700,
		com = 	{ 0.021735, 0.000000, -0.050000,},
		inertia = 
			{{ 0.003551, 0.000000, 0.000000,},
			{ 0.000000, 0.003092, 0.000000,},
			{ 0.000000, 0.000000, 0.000767,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.154744, 0.066319, -0.033159,},
		R_FM2 = 	{ 0.221063, 0.000000, -0.033159,},
		R_FM5 = 	{ 0.154744, -0.066319, -0.033159,},
		R_FCC = 	{ -0.066319, 0.000000, -0.033159,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.221063, 0.117000, 0.100000,},
		mesh_center = 	{ 0.043531, 0.000000, -0.050000,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.089000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 9.311400,
		com = 	{ 0.000000, 0.000000, -0.128883,},
		inertia = 
			{{ 0.161422, 0.000000, 0.000000,},
			{ 0.000000, 0.157077, 0.000000,},
			{ 0.000000, 0.000000, 0.031113,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.071364, -0.117750,},
		L_FLE = 	{ 0.000000, 0.053523, -0.356818,},
		L_FME = 	{ 0.000000, -0.053523, -0.356818,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.142727, 0.142727, 0.428182,},
		mesh_center = 	{ 0.000000, 0.000000, -0.178409,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.356818,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 3.030300,
		com = 	{ 0.000000, 0.000000, -0.184828,},
		inertia = 
			{{ 0.038964, 0.000000, 0.000000,},
			{ 0.000000, 0.037805, 0.000000,},
			{ 0.000000, 0.000000, 0.004626,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.084939, -0.021235,},
		L_TTC = 	{ 0.021235, 0.000000, -0.063704,},
		L_FAL = 	{ 0.000000, 0.063704, -0.424696,},
		L_TAM = 	{ 0.000000, -0.063704, -0.424696,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.127409, 0.127409, 0.509636,},
		mesh_center = 	{ 0.000000, 0.000000, -0.212348,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.424696,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.812700,
		com = 	{ 0.021735, 0.000000, -0.050000,},
		inertia = 
			{{ 0.003551, 0.000000, 0.000000,},
			{ 0.000000, 0.003092, 0.000000,},
			{ 0.000000, 0.000000, 0.000767,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.154744, -0.066319, -0.033159,},
		L_FM2 = 	{ 0.221063, 0.000000, -0.033159,},
		L_FM5 = 	{ 0.154744, 0.066319, -0.033159,},
		L_FCC = 	{ -0.066319, 0.000000, -0.033159,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.221063, 0.117000, 0.100000,},
		mesh_center = 	{ 0.043531, 0.000000, -0.050000,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.175746,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 9.229500,
		com = 	{ 0.000000, 0.000000, 0.109097,},
		inertia = 
			{{ 0.068384, 0.000000, 0.000000,},
			{ 0.000000, 0.045707, 0.000000,},
			{ 0.000000, 0.000000, 0.062816,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.116000, 0.000000, 0.099396,},
		MAI = 	{ -0.116000, 0.000000, 0.099396,},
		LV1 = 	{ -0.116000, 0.000000, 0.000000,},
		LV3 = 	{ -0.116000, 0.000000, -0.049698,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.232000, 0.234000, 0.298188,},
		mesh_center = 	{ 0.000000, 0.000000, 0.049698,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.198792,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 9.733500,
		com = 	{ 0.000000, 0.000000, 0.109282,},
		inertia = 
			{{ 0.103022, 0.000000, 0.000000,},
			{ 0.000000, 0.046775, 0.000000,},
			{ 0.000000, 0.000000, 0.095643,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.092800, 0.000000, 0.176618,},
		SJN = 	{ 0.097440, 0.000000, 0.132463,},
		TV2 = 	{ -0.092800, 0.000000, 0.132463,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.232000, 0.243600, 0.231811,},
		mesh_center = 	{ 0.000000, 0.000000, 0.088309,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.220772,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 4.208400,
		com = 	{ 0.000000, 0.000000, 0.121739,},
		inertia = 
			{{ 0.017210, 0.000000, 0.000000,},
			{ 0.000000, 0.020393, 0.000000,},
			{ 0.000000, 0.000000, 0.015964,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.070792, 0.141585,},
		R_HEAD = 	{ 0.000000, -0.070792, 0.141585,},
		SGL    = 	{ 0.070792, 0.000000, 0.141585,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.188780, 0.188780, 0.247773,},
		mesh_center = 	{ 0.000000, 0.000000, 0.094390,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.100500, 0.189772,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.606500,
		com = 	{ 0.000000, 0.000000, -0.153275,},
		inertia = 
			{{ 0.008810, 0.000000, 0.000000,},
			{ 0.000000, 0.007706, 0.000000,},
			{ 0.000000, 0.000000, 0.002497,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.026638,},
		R_HUM = 	{ 0.000000, -0.039957, -0.178474,},
		R_HLE = 	{ 0.000000, -0.039957, -0.266379,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.133190, 0.106552, 0.293017,},
		mesh_center = 	{ 0.000000, 0.000000, -0.133190,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.266379,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.869400,
		com = 	{ 0.000000, 0.000000, -0.116675,},
		inertia = 
			{{ 0.003879, 0.000000, 0.000000,},
			{ 0.000000, 0.003761, 0.000000,},
			{ 0.000000, 0.000000, 0.000503,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ -0.051184, 0.000000, -0.255922,},
		R_RSP = 	{ 0.051184, 0.000000, -0.255922,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.102369, 0.076776, 0.281514,},
		mesh_center = 	{ 0.000000, 0.000000, -0.127961,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.255922,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.352800,
		com = 	{ 0.000000, 0.000000, -0.056445,},
		inertia = 
			{{ 0.000570, 0.000000, 0.000000,},
			{ 0.000000, 0.000414, 0.000000,},
			{ 0.000000, 0.000000, 0.000227,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.131766,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.115295, 0.049412, 0.164708,},
		mesh_center = 	{ 0.000000, 0.000000, -0.082354,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.100500, 0.189772,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.606500,
		com = 	{ 0.000000, 0.000000, -0.153275,},
		inertia = 
			{{ 0.008810, 0.000000, 0.000000,},
			{ 0.000000, 0.007706, 0.000000,},
			{ 0.000000, 0.000000, 0.002497,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.026638,},
		L_HUM = 	{ 0.000000, 0.039957, -0.133190,},
		L_HLE = 	{ 0.000000, 0.039957, -0.266379,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.133190, 0.106552, 0.293017,},
		mesh_center = 	{ 0.000000, -0.000000, -0.133190,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.266379,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.869400,
		com = 	{ 0.000000, 0.000000, -0.116675,},
		inertia = 
			{{ 0.003879, 0.000000, 0.000000,},
			{ 0.000000, 0.003761, 0.000000,},
			{ 0.000000, 0.000000, 0.000503,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.051184, 0.000000, -0.255922,},
		L_RSP = 	{ 0.051184, 0.000000, -0.255922,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.102369, 0.076776, 0.281514,},
		mesh_center = 	{ 0.000000, 0.000000, -0.127961,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.255922,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.352800,
		com = 	{ 0.000000, 0.000000, -0.056445,},
		inertia = 
			{{ 0.000570, 0.000000, 0.000000,},
			{ 0.000000, 0.000414, 0.000000,},
			{ 0.000000, 0.000000, 0.000227,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.131766,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.115295, 0.049412, 0.164708,},
		mesh_center = 	{ 0.000000, 0.000000, -0.082354,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}