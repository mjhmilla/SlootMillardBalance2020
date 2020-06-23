return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {65.0},
	 subject_height = {1.63},
	 subject_weight = {96.00},
	 subject_gender = {"female"},
	 subject_pelvisWidth = {0.3150},
	 subject_hipCenterWidth = {0.2390},
	 subject_shoulderCenterWidth = {0.2430},
	 subject_heelAnkleXOffset = {0.0810},
	 subject_heelAnkleZOffset = {0.0850},
	 subject_shoulderNeckZOffset = {0.0360},
	 subject_footWidth = {0.1160},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Thigh_Proximal_R", body = "Thigh_R", point = {-0.069240, 0.000000, 0.034620,},},
	{name = "Thigh_Distal_R", body = "Thigh_R", point = {-0.069240, 0.000000, -0.086550,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.081000, 0.040600, -0.085000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.081000, -0.040600, -0.085000,},},
	{name = "ForeFoot_Medial_R", body = "Foot_R", point = {0.165656, 0.052200, -0.085000,},},
	{name = "ForeFoot_Lateral_R", body = "Foot_R", point = {0.165656, -0.052200, -0.085000,},},
	{name = "Thigh_Proximal_L", body = "Thigh_L", point = {-0.069240, 0.000000, 0.034620,},},
	{name = "Thigh_Distal_L", body = "Thigh_L", point = {-0.069240, 0.000000, -0.086550,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.081000, -0.040600, -0.085000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.081000, 0.040600, -0.085000,},},
	{name = "ForeFoot_Medial_L", body = "Foot_L", point = {0.165656, -0.052200, -0.085000,},},
	{name = "ForeFoot_Lateral_L", body = "Foot_L", point = {0.165656, 0.052200, -0.085000,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.031961, 0.007990, -0.079903,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.031961, 0.007990, -0.079903,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.031961, -0.007990, -0.079903,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.031961, -0.007990, -0.079903,},},
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
		mass   = 11.971200,
		com = 	{ 0.000000, 0.000000, 0.086622,},
		inertia = 
			{{ 0.065259, 0.000000, 0.000000,},
			{ 0.000000, 0.056250, 0.000000,},
			{ 0.000000, 0.000000, 0.068617,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.127887, 0.085258, 0.102310,},
		R_IAS = 	{ 0.127887, -0.085258, 0.102310,},
		L_IPS = 	{ -0.119361, 0.051155, 0.119361,},
		R_IPS = 	{ -0.119361, -0.051155, 0.119361,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.283500, 0.315000, 0.255774,},
		mesh_center = 	{ 0.000000, 0.000000, 0.042629,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.119500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 14.188800,
		com = 	{ 0.000000, 0.000000, -0.125047,},
		inertia = 
			{{ 0.231553, 0.000000, 0.000000,},
			{ 0.000000, 0.225320, 0.000000,},
			{ 0.000000, 0.000000, 0.044630,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.069240, -0.114246,},
		R_FLE = 	{ 0.000000, -0.051930, -0.346199,},
		R_FME = 	{ 0.000000, 0.051930, -0.346199,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.138480, 0.138480, 0.415439,},
		mesh_center = 	{ 0.000000, 0.000000, -0.173099,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.346199,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 4.617600,
		com = 	{ 0.000000, 0.000000, -0.179327,},
		inertia = 
			{{ 0.055892, 0.000000, 0.000000,},
			{ 0.000000, 0.054230, 0.000000,},
			{ 0.000000, 0.000000, 0.006636,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.082411, -0.020603,},
		R_TTC = 	{ 0.020603, 0.000000, -0.061808,},
		R_FAL = 	{ 0.000000, -0.061808, -0.412056,},
		R_TAM = 	{ 0.000000, 0.061808, -0.412056,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.123617, 0.123617, 0.494468,},
		mesh_center = 	{ 0.000000, 0.000000, -0.206028,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.412056,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.238400,
		com = 	{ 0.005094, 0.000000, -0.042500,},
		inertia = 
			{{ 0.005093, 0.000000, 0.000000,},
			{ 0.000000, 0.004435, 0.000000,},
			{ 0.000000, 0.000000, 0.001101,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.150139, 0.064345, -0.064345,},
		R_FM2 = 	{ 0.214484, 0.000000, -0.053621,},
		R_FM5 = 	{ 0.150139, -0.064345, -0.064345,},
		R_FCC = 	{ -0.064345, 0.000000, -0.021448,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.278829, 0.116000, 0.063750,},
		mesh_center = 	{ 0.058414, 0.000000, -0.053125,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.119500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 14.188800,
		com = 	{ 0.000000, 0.000000, -0.125047,},
		inertia = 
			{{ 0.231553, 0.000000, 0.000000,},
			{ 0.000000, 0.225320, 0.000000,},
			{ 0.000000, 0.000000, 0.044630,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.069240, -0.114246,},
		L_FLE = 	{ 0.000000, 0.051930, -0.346199,},
		L_FME = 	{ 0.000000, -0.051930, -0.346199,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.138480, 0.138480, 0.415439,},
		mesh_center = 	{ 0.000000, 0.000000, -0.173099,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.346199,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 4.617600,
		com = 	{ 0.000000, 0.000000, -0.179327,},
		inertia = 
			{{ 0.055892, 0.000000, 0.000000,},
			{ 0.000000, 0.054230, 0.000000,},
			{ 0.000000, 0.000000, 0.006636,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.082411, -0.020603,},
		L_TTC = 	{ 0.020603, 0.000000, -0.061808,},
		L_FAL = 	{ 0.000000, 0.061808, -0.412056,},
		L_TAM = 	{ 0.000000, -0.061808, -0.412056,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.123617, 0.123617, 0.494468,},
		mesh_center = 	{ 0.000000, 0.000000, -0.206028,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.412056,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.238400,
		com = 	{ 0.005094, 0.000000, -0.042500,},
		inertia = 
			{{ 0.005093, 0.000000, 0.000000,},
			{ 0.000000, 0.004435, 0.000000,},
			{ 0.000000, 0.000000, 0.001101,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.150139, -0.064345, -0.064345,},
		L_FM2 = 	{ 0.214484, 0.000000, -0.053621,},
		L_FM5 = 	{ 0.150139, 0.064345, -0.064345,},
		L_FCC = 	{ -0.064345, 0.000000, -0.021448,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.278829, 0.116000, 0.063750,},
		mesh_center = 	{ 0.058414, 0.000000, -0.053125,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.170516,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 14.064000,
		com = 	{ 0.000000, 0.000000, 0.105850,},
		inertia = 
			{{ 0.098093, 0.000000, 0.000000,},
			{ 0.000000, 0.065565, 0.000000,},
			{ 0.000000, 0.000000, 0.090107,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.116000, 0.000000, 0.096438,},
		MAI = 	{ -0.116000, 0.000000, 0.096438,},
		LV1 = 	{ -0.116000, 0.000000, 0.000000,},
		LV3 = 	{ -0.116000, 0.000000, -0.048219,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.232000, 0.315000, 0.289313,},
		mesh_center = 	{ 0.000000, 0.000000, 0.048219,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.192876,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 14.832000,
		com = 	{ 0.000000, 0.000000, 0.106030,},
		inertia = 
			{{ 0.147781, 0.000000, 0.000000,},
			{ 0.000000, 0.067097, 0.000000,},
			{ 0.000000, 0.000000, 0.137195,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.092800, 0.000000, 0.171361,},
		SJN = 	{ 0.097440, 0.000000, 0.128521,},
		TV2 = 	{ -0.092800, 0.000000, 0.128521,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.232000, 0.255150, 0.224912,},
		mesh_center = 	{ 0.000000, 0.000000, 0.085681,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.214202,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 6.412800,
		com = 	{ 0.000000, 0.000000, 0.118116,},
		inertia = 
			{{ 0.024687, 0.000000, 0.000000,},
			{ 0.000000, 0.029254, 0.000000,},
			{ 0.000000, 0.000000, 0.022899,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.068685, 0.137371,},
		R_HEAD = 	{ 0.000000, -0.068685, 0.137371,},
		SGL    = 	{ 0.068685, 0.000000, 0.137371,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.183161, 0.183161, 0.240399,},
		mesh_center = 	{ 0.000000, 0.000000, 0.091581,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.121500, 0.178202,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.448000,
		com = 	{ 0.000000, 0.000000, -0.148713,},
		inertia = 
			{{ 0.012637, 0.000000, 0.000000,},
			{ 0.000000, 0.011054, 0.000000,},
			{ 0.000000, 0.000000, 0.003582,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.025845,},
		R_HUM = 	{ 0.000000, -0.038768, -0.173162,},
		R_HLE = 	{ 0.000000, -0.038768, -0.258451,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.129226, 0.103381, 0.284296,},
		mesh_center = 	{ 0.000000, 0.000000, -0.129226,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.258451,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.324800,
		com = 	{ 0.000000, 0.000000, -0.113202,},
		inertia = 
			{{ 0.005564, 0.000000, 0.000000,},
			{ 0.000000, 0.005395, 0.000000,},
			{ 0.000000, 0.000000, 0.000722,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ -0.049661, 0.000000, -0.248305,},
		R_RSP = 	{ 0.049661, 0.000000, -0.248305,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.099322, 0.074491, 0.273135,},
		mesh_center = 	{ 0.000000, 0.000000, -0.124152,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.248305,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.537600,
		com = 	{ 0.000000, 0.000000, -0.054765,},
		inertia = 
			{{ 0.000817, 0.000000, 0.000000,},
			{ 0.000000, 0.000594, 0.000000,},
			{ 0.000000, 0.000000, 0.000326,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.127845,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.111864, 0.047942, 0.159806,},
		mesh_center = 	{ 0.000000, 0.000000, -0.079903,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.121500, 0.178202,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.448000,
		com = 	{ 0.000000, 0.000000, -0.148713,},
		inertia = 
			{{ 0.012637, 0.000000, 0.000000,},
			{ 0.000000, 0.011054, 0.000000,},
			{ 0.000000, 0.000000, 0.003582,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.025845,},
		L_HUM = 	{ 0.000000, 0.038768, -0.129226,},
		L_HLE = 	{ 0.000000, 0.038768, -0.258451,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.129226, 0.103381, 0.284296,},
		mesh_center = 	{ 0.000000, -0.000000, -0.129226,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.258451,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.324800,
		com = 	{ 0.000000, 0.000000, -0.113202,},
		inertia = 
			{{ 0.005564, 0.000000, 0.000000,},
			{ 0.000000, 0.005395, 0.000000,},
			{ 0.000000, 0.000000, 0.000722,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.049661, 0.000000, -0.248305,},
		L_RSP = 	{ 0.049661, 0.000000, -0.248305,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.099322, 0.074491, 0.273135,},
		mesh_center = 	{ 0.000000, 0.000000, -0.124152,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.248305,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.537600,
		com = 	{ 0.000000, 0.000000, -0.054765,},
		inertia = 
			{{ 0.000817, 0.000000, 0.000000,},
			{ 0.000000, 0.000594, 0.000000,},
			{ 0.000000, 0.000000, 0.000326,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.127845,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.111864, 0.047942, 0.159806,},
		mesh_center = 	{ 0.000000, 0.000000, -0.079903,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}