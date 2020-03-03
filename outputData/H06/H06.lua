return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {25.0},
	 subject_height = {1.50},
	 subject_weight = {55.00},
	 subject_gender = {"female"},
	 subject_pelvisWidth = {0.1810},
	 subject_hipCenterWidth = {0.1370},
	 subject_shoulderCenterWidth = {0.2080},
	 subject_heelAnkleXOffset = {0.0950},
	 subject_heelAnkleZOffset = {0.1000},
	 subject_shoulderNeckZOffset = {0.0420},
	 subject_footWidth = {0.1140},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Pelvis_L", body = "Pelvis", point = {0.000000, 0.156916, 0.109841,},},
	{name = "Pelvis_R", body = "Pelvis", point = {0.000000, -0.156916, 0.109841,},},
	{name = "Pelvis_Back", body = "Pelvis", point = {-0.117687, 0.000000, 0.109841,},},
	{name = "Pelvis_Front", body = "Pelvis", point = {0.117687, 0.000000, 0.109841,},},
	{name = "Thigh_R", body = "Thigh_R", point = {0.063718, -0.000000, -0.238941,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.095000, 0.057000, -0.100000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.095000, -0.057000, -0.100000,},},
	{name = "Toe_R", body = "Foot_R", point = {0.053033, 0.000000, -0.100000,},},
	{name = "Thigh_L", body = "Thigh_L", point = {0.063718, 0.000000, -0.238941,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.095000, -0.057000, -0.100000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.095000, 0.057000, -0.100000,},},
	{name = "Toe_L", body = "Foot_L", point = {0.053033, 0.000000, -0.100000,},},
	{name = "UpperTrunk_Front", body = "UpperTrunk", point = {0.078847, 0.000000, 0.059135,},},
	{name = "UpperTrunk_Back", body = "UpperTrunk", point = {-0.088703, 0.000000, 0.078847,},},
	{name = "ProximalMetacarpal_Medial_R", body = "Hand_R", point = {-0.029412, 0.022059, -0.029412,},},
	{name = "ProximalMetacarpal_Lateral_R", body = "Hand_R", point = {0.029412, 0.022059, -0.029412,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.029412, 0.022059, -0.088236,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.029412, 0.022059, -0.088236,},},
	{name = "ProximalMetacarpal_Medial_L", body = "Hand_L", point = {-0.029412, -0.022059, -0.029412,},},
	{name = "ProximalMetacarpal_Lateral_L", body = "Hand_L", point = {0.029412, -0.022059, -0.029412,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.029412, -0.022059, -0.088236,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.029412, -0.022059, -0.088236,},},
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
		mass   = 6.858500,
		com = 	{ 0.000000, 0.000000, 0.079714,},
		inertia = 
			{{ 0.031662, 0.000000, 0.000000,},
			{ 0.000000, 0.027291, 0.000000,},
			{ 0.000000, 0.000000, 0.033291,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.117687, 0.078458, 0.094150,},
		R_IAS = 	{ 0.117687, -0.078458, 0.094150,},
		L_IPS = 	{ -0.109841, 0.109841, 0.109841,},
		R_IPS = 	{ -0.109841, -0.109841, 0.109841,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.181000, 0.226250, 0.235375,},
		mesh_center = 	{ 0.000000, 0.000000, 0.039229,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.068500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 8.129000,
		com = 	{ 0.000000, 0.000000, -0.115074,},
		inertia = 
			{{ 0.112344, 0.000000, 0.000000,},
			{ 0.000000, 0.109320, 0.000000,},
			{ 0.000000, 0.000000, 0.021653,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.063718, -0.105134,},
		R_FLE = 	{ 0.000000, -0.047788, -0.318588,},
		R_FME = 	{ 0.000000, 0.047788, -0.318588,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.159294, 0.127435, 0.382305,},
		mesh_center = 	{ 0.000000, 0.000000, -0.159294,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.318588,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.645500,
		com = 	{ 0.000000, 0.000000, -0.165025,},
		inertia = 
			{{ 0.027118, 0.000000, 0.000000,},
			{ 0.000000, 0.026311, 0.000000,},
			{ 0.000000, 0.000000, 0.003220,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.075839, -0.018960,},
		R_TTC = 	{ 0.018960, 0.000000, -0.056879,},
		R_FAL = 	{ 0.000000, -0.056879, -0.379193,},
		R_TAM = 	{ 0.000000, 0.056879, -0.379193,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.113758, 0.113758, 0.455032,},
		mesh_center = 	{ 0.000000, 0.000000, -0.189597,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.379193,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.709500,
		com = 	{ -0.015773, 0.000000, -0.050000,},
		inertia = 
			{{ 0.002471, 0.000000, 0.000000,},
			{ 0.000000, 0.002152, 0.000000,},
			{ 0.000000, 0.000000, 0.000534,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.138164, 0.019738, -0.029607,},
		R_FM2 = 	{ 0.138164, 0.000000, -0.029607,},
		R_FM5 = 	{ 0.138164, -0.019738, -0.029607,},
		R_FCC = 	{ -0.059213, 0.000000, -0.029607,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.197378, 0.114000, 0.100000,},
		mesh_center = 	{ 0.003689, 0.000000, -0.050000,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.068500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 8.129000,
		com = 	{ 0.000000, 0.000000, -0.115074,},
		inertia = 
			{{ 0.112344, 0.000000, 0.000000,},
			{ 0.000000, 0.109320, 0.000000,},
			{ 0.000000, 0.000000, 0.021653,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.063718, -0.105134,},
		L_FLE = 	{ 0.000000, 0.047788, -0.318588,},
		L_FME = 	{ 0.000000, -0.047788, -0.318588,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.159294, 0.127435, 0.382305,},
		mesh_center = 	{ 0.000000, 0.000000, -0.159294,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.318588,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.645500,
		com = 	{ 0.000000, 0.000000, -0.165025,},
		inertia = 
			{{ 0.027118, 0.000000, 0.000000,},
			{ 0.000000, 0.026311, 0.000000,},
			{ 0.000000, 0.000000, 0.003220,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.075839, -0.018960,},
		L_TTC = 	{ 0.018960, 0.000000, -0.056879,},
		L_FAL = 	{ 0.000000, 0.056879, -0.379193,},
		L_TAM = 	{ 0.000000, -0.056879, -0.379193,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.113758, 0.113758, 0.455032,},
		mesh_center = 	{ 0.000000, 0.000000, -0.189597,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.379193,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.709500,
		com = 	{ -0.015773, 0.000000, -0.050000,},
		inertia = 
			{{ 0.002471, 0.000000, 0.000000,},
			{ 0.000000, 0.002152, 0.000000,},
			{ 0.000000, 0.000000, 0.000534,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.138164, -0.019738, -0.029607,},
		L_FM2 = 	{ 0.138164, 0.000000, -0.029607,},
		L_FM5 = 	{ 0.138164, 0.019738, -0.029607,},
		L_FCC = 	{ -0.059213, 0.000000, -0.029607,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.197378, 0.114000, 0.100000,},
		mesh_center = 	{ 0.003689, 0.000000, -0.050000,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.156916,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 8.057500,
		com = 	{ 0.000000, 0.000000, 0.097408,},
		inertia = 
			{{ 0.047592, 0.000000, 0.000000,},
			{ 0.000000, 0.031810, 0.000000,},
			{ 0.000000, 0.000000, 0.043718,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.088746, 0.000000, 0.088746,},
		MAI = 	{ -0.088746, 0.000000, 0.088746,},
		LV1 = 	{ -0.088746, 0.000000, 0.000000,},
		LV3 = 	{ -0.088746, 0.000000, -0.044373,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.162900, 0.199100, 0.266239,},
		mesh_center = 	{ 0.000000, 0.000000, 0.044373,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.177493,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 8.497500,
		com = 	{ 0.000000, 0.000000, 0.097573,},
		inertia = 
			{{ 0.071700, 0.000000, 0.000000,},
			{ 0.000000, 0.032554, 0.000000,},
			{ 0.000000, 0.000000, 0.066564,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.078847, 0.000000, 0.157695,},
		SJN = 	{ 0.082790, 0.000000, 0.118271,},
		TV2 = 	{ -0.078847, 0.000000, 0.118271,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.136500, 0.260000, 0.206974,},
		mesh_center = 	{ 0.000000, 0.000000, 0.078847,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.197118,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 3.674000,
		com = 	{ 0.000000, 0.000000, 0.108696,},
		inertia = 
			{{ 0.011978, 0.000000, 0.000000,},
			{ 0.000000, 0.014193, 0.000000,},
			{ 0.000000, 0.000000, 0.011110,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.063207, 0.126415,},
		R_HEAD = 	{ 0.000000, -0.063207, 0.126415,},
		SGL    = 	{ 0.063207, 0.000000, 0.126415,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.168553, 0.168553, 0.221226,},
		mesh_center = 	{ 0.000000, 0.000000, 0.084277,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.104000, 0.155118,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.402500,
		com = 	{ 0.000000, 0.000000, -0.136852,},
		inertia = 
			{{ 0.006131, 0.000000, 0.000000,},
			{ 0.000000, 0.005363, 0.000000,},
			{ 0.000000, 0.000000, 0.001738,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.023784,},
		R_HUM = 	{ 0.000000, -0.035676, -0.159352,},
		R_HLE = 	{ 0.000000, -0.035676, -0.237839,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.118919, 0.095135, 0.261622,},
		mesh_center = 	{ 0.000000, 0.000000, -0.118919,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.237839,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.759000,
		com = 	{ 0.000000, 0.000000, -0.104174,},
		inertia = 
			{{ 0.002700, 0.000000, 0.000000,},
			{ 0.000000, 0.002617, 0.000000,},
			{ 0.000000, 0.000000, 0.000350,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ 0.045700, 0.000000, -0.228501,},
		R_RSP = 	{ -0.045700, 0.000000, -0.228501,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.091401, 0.068550, 0.251352,},
		mesh_center = 	{ 0.000000, 0.000000, -0.114251,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.228501,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.308000,
		com = 	{ 0.000000, 0.000000, -0.050398,},
		inertia = 
			{{ 0.000397, 0.000000, 0.000000,},
			{ 0.000000, 0.000288, 0.000000,},
			{ 0.000000, 0.000000, 0.000158,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.117648,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.102942, 0.044118, 0.147061,},
		mesh_center = 	{ 0.000000, 0.000000, -0.073530,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.104000, 0.155118,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.402500,
		com = 	{ 0.000000, 0.000000, -0.136852,},
		inertia = 
			{{ 0.006131, 0.000000, 0.000000,},
			{ 0.000000, 0.005363, 0.000000,},
			{ 0.000000, 0.000000, 0.001738,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.023784,},
		L_HUM = 	{ 0.000000, 0.035676, -0.118919,},
		L_HLE = 	{ 0.000000, 0.035676, -0.237839,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.118919, 0.095135, 0.261622,},
		mesh_center = 	{ 0.000000, -0.000000, -0.118919,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.237839,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.759000,
		com = 	{ 0.000000, 0.000000, -0.104174,},
		inertia = 
			{{ 0.002700, 0.000000, 0.000000,},
			{ 0.000000, 0.002617, 0.000000,},
			{ 0.000000, 0.000000, 0.000350,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.045700, 0.000000, -0.228501,},
		L_RSP = 	{ 0.045700, 0.000000, -0.228501,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.091401, 0.068550, 0.251352,},
		mesh_center = 	{ 0.000000, 0.000000, -0.114251,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.228501,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.308000,
		com = 	{ 0.000000, 0.000000, -0.050398,},
		inertia = 
			{{ 0.000397, 0.000000, 0.000000,},
			{ 0.000000, 0.000288, 0.000000,},
			{ 0.000000, 0.000000, 0.000158,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.117648,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.102942, 0.044118, 0.147061,},
		mesh_center = 	{ 0.000000, 0.000000, -0.073530,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}