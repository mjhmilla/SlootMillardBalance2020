return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {25.0},
	 subject_height = {1.61},
	 subject_weight = {52.00},
	 subject_gender = {"female"},
	 subject_pelvisWidth = {0.2530},
	 subject_hipCenterWidth = {0.1920},
	 subject_shoulderCenterWidth = {0.2260},
	 subject_heelAnkleXOffset = {0.1020},
	 subject_heelAnkleZOffset = {0.0980},
	 subject_shoulderNeckZOffset = {0.0420},
	 subject_footWidth = {0.1180},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Pelvis_L", body = "Pelvis", point = {0.000000, 0.168424, 0.117897,},},
	{name = "Pelvis_R", body = "Pelvis", point = {0.000000, -0.168424, 0.117897,},},
	{name = "Pelvis_Back", body = "Pelvis", point = {-0.126318, 0.000000, 0.117897,},},
	{name = "Pelvis_Front", body = "Pelvis", point = {0.126318, 0.000000, 0.117897,},},
	{name = "Thigh_R", body = "Thigh_R", point = {0.068390, -0.000000, -0.256463,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.102000, 0.059000, -0.098000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.102000, -0.059000, -0.098000,},},
	{name = "Toe_R", body = "Foot_R", point = {0.056889, 0.000000, -0.098000,},},
	{name = "Thigh_L", body = "Thigh_L", point = {0.068390, 0.000000, -0.256463,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.102000, -0.059000, -0.098000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.102000, 0.059000, -0.098000,},},
	{name = "Toe_L", body = "Foot_L", point = {0.056889, 0.000000, -0.098000,},},
	{name = "UpperTrunk_Front", body = "UpperTrunk", point = {0.084629, 0.000000, 0.063472,},},
	{name = "UpperTrunk_Back", body = "UpperTrunk", point = {-0.095208, 0.000000, 0.084629,},},
	{name = "ProximalMetacarpal_Medial_R", body = "Hand_R", point = {-0.031569, 0.023677, -0.031569,},},
	{name = "ProximalMetacarpal_Lateral_R", body = "Hand_R", point = {0.031569, 0.023677, -0.031569,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.031569, 0.023677, -0.094707,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.031569, 0.023677, -0.094707,},},
	{name = "ProximalMetacarpal_Medial_L", body = "Hand_L", point = {-0.031569, -0.023677, -0.031569,},},
	{name = "ProximalMetacarpal_Lateral_L", body = "Hand_L", point = {0.031569, -0.023677, -0.031569,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.031569, -0.023677, -0.094707,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.031569, -0.023677, -0.094707,},},
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
		mass   = 6.484400,
		com = 	{ 0.000000, 0.000000, 0.085559,},
		inertia = 
			{{ 0.034487, 0.000000, 0.000000,},
			{ 0.000000, 0.029725, 0.000000,},
			{ 0.000000, 0.000000, 0.036261,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.126318, 0.084212, 0.101054,},
		R_IAS = 	{ 0.126318, -0.084212, 0.101054,},
		L_IPS = 	{ -0.117897, 0.117897, 0.117897,},
		R_IPS = 	{ -0.117897, -0.117897, 0.117897,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.253000, 0.316250, 0.252635,},
		mesh_center = 	{ 0.000000, 0.000000, 0.042106,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.096000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 7.685600,
		com = 	{ 0.000000, 0.000000, -0.123513,},
		inertia = 
			{{ 0.122365, 0.000000, 0.000000,},
			{ 0.000000, 0.119072, 0.000000,},
			{ 0.000000, 0.000000, 0.023585,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.068390, -0.112844,},
		R_FLE = 	{ 0.000000, -0.051293, -0.341951,},
		R_FME = 	{ 0.000000, 0.051293, -0.341951,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.170976, 0.136780, 0.410341,},
		mesh_center = 	{ 0.000000, 0.000000, -0.170976,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.341951,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.501200,
		com = 	{ 0.000000, 0.000000, -0.177127,},
		inertia = 
			{{ 0.029537, 0.000000, 0.000000,},
			{ 0.000000, 0.028658, 0.000000,},
			{ 0.000000, 0.000000, 0.003507,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.081400, -0.020350,},
		R_TTC = 	{ 0.020350, 0.000000, -0.061050,},
		R_FAL = 	{ 0.000000, -0.061050, -0.407001,},
		R_TAM = 	{ 0.000000, 0.061050, -0.407001,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.122100, 0.122100, 0.488401,},
		mesh_center = 	{ 0.000000, 0.000000, -0.203500,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.407001,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.670800,
		com = 	{ -0.016963, 0.000000, -0.049000,},
		inertia = 
			{{ 0.002692, 0.000000, 0.000000,},
			{ 0.000000, 0.002344, 0.000000,},
			{ 0.000000, 0.000000, 0.000582,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.148296, 0.021185, -0.031778,},
		R_FM2 = 	{ 0.148296, 0.000000, -0.031778,},
		R_FM5 = 	{ 0.148296, -0.021185, -0.031778,},
		R_FCC = 	{ -0.063556, 0.000000, -0.031778,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.211852, 0.118000, 0.098000,},
		mesh_center = 	{ 0.003926, 0.000000, -0.049000,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.096000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 7.685600,
		com = 	{ 0.000000, 0.000000, -0.123513,},
		inertia = 
			{{ 0.122365, 0.000000, 0.000000,},
			{ 0.000000, 0.119072, 0.000000,},
			{ 0.000000, 0.000000, 0.023585,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.068390, -0.112844,},
		L_FLE = 	{ 0.000000, 0.051293, -0.341951,},
		L_FME = 	{ 0.000000, -0.051293, -0.341951,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.170976, 0.136780, 0.410341,},
		mesh_center = 	{ 0.000000, 0.000000, -0.170976,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.341951,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.501200,
		com = 	{ 0.000000, 0.000000, -0.177127,},
		inertia = 
			{{ 0.029537, 0.000000, 0.000000,},
			{ 0.000000, 0.028658, 0.000000,},
			{ 0.000000, 0.000000, 0.003507,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.081400, -0.020350,},
		L_TTC = 	{ 0.020350, 0.000000, -0.061050,},
		L_FAL = 	{ 0.000000, 0.061050, -0.407001,},
		L_TAM = 	{ 0.000000, -0.061050, -0.407001,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.122100, 0.122100, 0.488401,},
		mesh_center = 	{ 0.000000, 0.000000, -0.203500,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.407001,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.670800,
		com = 	{ -0.016963, 0.000000, -0.049000,},
		inertia = 
			{{ 0.002692, 0.000000, 0.000000,},
			{ 0.000000, 0.002344, 0.000000,},
			{ 0.000000, 0.000000, 0.000582,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.148296, -0.021185, -0.031778,},
		L_FM2 = 	{ 0.148296, 0.000000, -0.031778,},
		L_FM5 = 	{ 0.148296, 0.021185, -0.031778,},
		L_FCC = 	{ -0.063556, 0.000000, -0.031778,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.211852, 0.118000, 0.098000,},
		mesh_center = 	{ 0.003926, 0.000000, -0.049000,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.168424,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 7.618000,
		com = 	{ 0.000000, 0.000000, 0.104551,},
		inertia = 
			{{ 0.051838, 0.000000, 0.000000,},
			{ 0.000000, 0.034648, 0.000000,},
			{ 0.000000, 0.000000, 0.047618,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.095254, 0.000000, 0.095254,},
		MAI = 	{ -0.095254, 0.000000, 0.095254,},
		LV1 = 	{ -0.095254, 0.000000, 0.000000,},
		LV3 = 	{ -0.095254, 0.000000, -0.047627,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.227700, 0.278300, 0.285763,},
		mesh_center = 	{ 0.000000, 0.000000, 0.047627,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.190509,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 8.034000,
		com = 	{ 0.000000, 0.000000, 0.104729,},
		inertia = 
			{{ 0.078096, 0.000000, 0.000000,},
			{ 0.000000, 0.035458, 0.000000,},
			{ 0.000000, 0.000000, 0.072502,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.084629, 0.000000, 0.169259,},
		SJN = 	{ 0.088861, 0.000000, 0.126944,},
		TV2 = 	{ -0.084629, 0.000000, 0.126944,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.148313, 0.282500, 0.222152,},
		mesh_center = 	{ 0.000000, 0.000000, 0.084629,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.211573,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 3.473600,
		com = 	{ 0.000000, 0.000000, 0.116667,},
		inertia = 
			{{ 0.013046, 0.000000, 0.000000,},
			{ 0.000000, 0.015459, 0.000000,},
			{ 0.000000, 0.000000, 0.012101,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.067843, 0.135685,},
		R_HEAD = 	{ 0.000000, -0.067843, 0.135685,},
		SGL    = 	{ 0.067843, 0.000000, 0.135685,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.180914, 0.180914, 0.237449,},
		mesh_center = 	{ 0.000000, 0.000000, 0.090457,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.113000, 0.169573,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.326000,
		com = 	{ 0.000000, 0.000000, -0.146888,},
		inertia = 
			{{ 0.006678, 0.000000, 0.000000,},
			{ 0.000000, 0.005841, 0.000000,},
			{ 0.000000, 0.000000, 0.001893,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.025528,},
		R_HUM = 	{ 0.000000, -0.038292, -0.171038,},
		R_HLE = 	{ 0.000000, -0.038292, -0.255280,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.127640, 0.102112, 0.280808,},
		mesh_center = 	{ 0.000000, 0.000000, -0.127640,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.255280,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.717600,
		com = 	{ 0.000000, 0.000000, -0.111813,},
		inertia = 
			{{ 0.002940, 0.000000, 0.000000,},
			{ 0.000000, 0.002851, 0.000000,},
			{ 0.000000, 0.000000, 0.000381,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ 0.049052, 0.000000, -0.245258,},
		R_RSP = 	{ -0.049052, 0.000000, -0.245258,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.098103, 0.073577, 0.269784,},
		mesh_center = 	{ 0.000000, 0.000000, -0.122629,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.245258,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.291200,
		com = 	{ 0.000000, 0.000000, -0.054093,},
		inertia = 
			{{ 0.000432, 0.000000, 0.000000,},
			{ 0.000000, 0.000314, 0.000000,},
			{ 0.000000, 0.000000, 0.000172,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.126276,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.110491, 0.047353, 0.157845,},
		mesh_center = 	{ 0.000000, 0.000000, -0.078922,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.113000, 0.169573,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.326000,
		com = 	{ 0.000000, 0.000000, -0.146888,},
		inertia = 
			{{ 0.006678, 0.000000, 0.000000,},
			{ 0.000000, 0.005841, 0.000000,},
			{ 0.000000, 0.000000, 0.001893,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.025528,},
		L_HUM = 	{ 0.000000, 0.038292, -0.127640,},
		L_HLE = 	{ 0.000000, 0.038292, -0.255280,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.127640, 0.102112, 0.280808,},
		mesh_center = 	{ 0.000000, -0.000000, -0.127640,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.255280,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.717600,
		com = 	{ 0.000000, 0.000000, -0.111813,},
		inertia = 
			{{ 0.002940, 0.000000, 0.000000,},
			{ 0.000000, 0.002851, 0.000000,},
			{ 0.000000, 0.000000, 0.000381,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.049052, 0.000000, -0.245258,},
		L_RSP = 	{ 0.049052, 0.000000, -0.245258,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.098103, 0.073577, 0.269784,},
		mesh_center = 	{ 0.000000, 0.000000, -0.122629,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.245258,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.291200,
		com = 	{ 0.000000, 0.000000, -0.054093,},
		inertia = 
			{{ 0.000432, 0.000000, 0.000000,},
			{ 0.000000, 0.000314, 0.000000,},
			{ 0.000000, 0.000000, 0.000172,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.126276,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.110491, 0.047353, 0.157845,},
		mesh_center = 	{ 0.000000, 0.000000, -0.078922,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}