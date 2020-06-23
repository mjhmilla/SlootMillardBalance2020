return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {65.0},
	 subject_height = {1.72},
	 subject_weight = {100.00},
	 subject_gender = {"male"},
	 subject_pelvisWidth = {0.2330},
	 subject_hipCenterWidth = {0.1770},
	 subject_shoulderCenterWidth = {0.2540},
	 subject_heelAnkleXOffset = {0.0920},
	 subject_heelAnkleZOffset = {0.1020},
	 subject_shoulderNeckZOffset = {0.0430},
	 subject_footWidth = {0.1260},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Thigh_Proximal_R", body = "Thigh_R", point = {-0.083421, 0.000000, 0.041711,},},
	{name = "Thigh_Distal_R", body = "Thigh_R", point = {-0.083421, 0.000000, -0.104277,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.092000, 0.044100, -0.102000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.092000, -0.044100, -0.102000,},},
	{name = "ForeFoot_Medial_R", body = "Foot_R", point = {0.201235, 0.056700, -0.102000,},},
	{name = "ForeFoot_Lateral_R", body = "Foot_R", point = {0.201235, -0.056700, -0.102000,},},
	{name = "Thigh_Proximal_L", body = "Thigh_L", point = {-0.083421, 0.000000, 0.041711,},},
	{name = "Thigh_Distal_L", body = "Thigh_L", point = {-0.083421, 0.000000, -0.104277,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.092000, -0.044100, -0.102000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.092000, 0.044100, -0.102000,},},
	{name = "ForeFoot_Medial_L", body = "Foot_L", point = {0.201235, -0.056700, -0.102000,},},
	{name = "ForeFoot_Lateral_L", body = "Foot_L", point = {0.201235, 0.056700, -0.102000,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.037127, 0.009282, -0.092817,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.037127, 0.009282, -0.092817,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.037127, -0.009282, -0.092817,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.037127, -0.009282, -0.092817,},},
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
		mass   = 11.170000,
		com = 	{ 0.000000, 0.000000, 0.055922,},
		inertia = 
			{{ 0.087535, 0.000000, 0.000000,},
			{ 0.000000, 0.070264, 0.000000,},
			{ 0.000000, 0.000000, 0.079746,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.107957, 0.071971, 0.086366,},
		R_IAS = 	{ 0.107957, -0.071971, 0.086366,},
		L_IPS = 	{ -0.100760, 0.043183, 0.100760,},
		R_IPS = 	{ -0.100760, -0.043183, 0.100760,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.209700, 0.233000, 0.215914,},
		mesh_center = 	{ 0.000000, 0.000000, 0.035986,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.088500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 14.160000,
		com = 	{ 0.000000, 0.000000, -0.170805,},
		inertia = 
			{{ 0.266656, 0.000000, 0.000000,},
			{ 0.000000, 0.266656, 0.000000,},
			{ 0.000000, 0.000000, 0.054693,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.083421, -0.137645,},
		R_FLE = 	{ 0.000000, -0.062566, -0.417107,},
		R_FME = 	{ 0.000000, 0.062566, -0.417107,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.166843, 0.166843, 0.500529,},
		mesh_center = 	{ 0.000000, 0.000000, -0.208554,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.417107,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 4.330000,
		com = 	{ 0.000000, 0.000000, -0.191178,},
		inertia = 
			{{ 0.051617, 0.000000, 0.000000,},
			{ 0.000000, 0.049581, 0.000000,},
			{ 0.000000, 0.000000, 0.008524,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.086998, -0.021749,},
		R_TTC = 	{ 0.021749, 0.000000, -0.065248,},
		R_FAL = 	{ 0.000000, -0.065248, -0.434989,},
		R_TAM = 	{ 0.000000, 0.065248, -0.434989,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.130497, 0.130497, 0.521987,},
		mesh_center = 	{ 0.000000, 0.000000, -0.217495,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.434989,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.370000,
		com = 	{ 0.020577, 0.000000, -0.051000,},
		inertia = 
			{{ 0.005883, 0.000000, 0.000000,},
			{ 0.000000, 0.005347, 0.000000,},
			{ 0.000000, 0.000000, 0.001370,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.178491, 0.076496, -0.076496,},
		R_FM2 = 	{ 0.254987, 0.000000, -0.063747,},
		R_FM5 = 	{ 0.178491, -0.076496, -0.076496,},
		R_FCC = 	{ -0.076496, 0.000000, -0.025499,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.331483, 0.126000, 0.076500,},
		mesh_center = 	{ 0.073741, 0.000000, -0.063750,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.088500, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 14.160000,
		com = 	{ 0.000000, 0.000000, -0.170805,},
		inertia = 
			{{ 0.266656, 0.000000, 0.000000,},
			{ 0.000000, 0.266656, 0.000000,},
			{ 0.000000, 0.000000, 0.054693,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.083421, -0.137645,},
		L_FLE = 	{ 0.000000, 0.062566, -0.417107,},
		L_FME = 	{ 0.000000, -0.062566, -0.417107,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.166843, 0.166843, 0.500529,},
		mesh_center = 	{ 0.000000, 0.000000, -0.208554,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.417107,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 4.330000,
		com = 	{ 0.000000, 0.000000, -0.191178,},
		inertia = 
			{{ 0.051617, 0.000000, 0.000000,},
			{ 0.000000, 0.049581, 0.000000,},
			{ 0.000000, 0.000000, 0.008524,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.086998, -0.021749,},
		L_TTC = 	{ 0.021749, 0.000000, -0.065248,},
		L_FAL = 	{ 0.000000, 0.065248, -0.434989,},
		L_TAM = 	{ 0.000000, -0.065248, -0.434989,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.130497, 0.130497, 0.521987,},
		mesh_center = 	{ 0.000000, 0.000000, -0.217495,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.434989,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.370000,
		com = 	{ 0.020577, 0.000000, -0.051000,},
		inertia = 
			{{ 0.005883, 0.000000, 0.000000,},
			{ 0.000000, 0.005347, 0.000000,},
			{ 0.000000, 0.000000, 0.001370,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.178491, -0.076496, -0.076496,},
		L_FM2 = 	{ 0.254987, 0.000000, -0.063747,},
		L_FM5 = 	{ 0.178491, 0.076496, -0.076496,},
		L_FCC = 	{ -0.076496, 0.000000, -0.025499,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.331483, 0.126000, 0.076500,},
		mesh_center = 	{ 0.073741, 0.000000, -0.063750,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.143943,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 16.330000,
		com = 	{ 0.000000, 0.000000, 0.117053,},
		inertia = 
			{{ 0.171963, 0.000000, 0.000000,},
			{ 0.000000, 0.108577, 0.000000,},
			{ 0.000000, 0.000000, 0.162118,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.151000, 0.000000, 0.106450,},
		MAI = 	{ -0.151000, 0.000000, 0.106450,},
		LV1 = 	{ -0.151000, 0.000000, 0.000000,},
		LV3 = 	{ -0.151000, 0.000000, -0.053225,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.302000, 0.302000, 0.319351,},
		mesh_center = 	{ 0.000000, 0.000000, 0.053225,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.212901,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 15.960000,
		com = 	{ 0.000000, 0.000000, 0.118011,},
		inertia = 
			{{ 0.232844, 0.000000, 0.000000,},
			{ 0.000000, 0.093494, 0.000000,},
			{ 0.000000, 0.000000, 0.197418,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.120800, 0.000000, 0.191344,},
		SJN = 	{ 0.126840, 0.000000, 0.143508,},
		TV2 = 	{ -0.120800, 0.000000, 0.143508,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.302000, 0.317100, 0.251139,},
		mesh_center = 	{ 0.000000, 0.000000, 0.095672,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.239180,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 6.940000,
		com = 	{ 0.000000, 0.000000, 0.119937,},
		inertia = 
			{{ 0.036691, 0.000000, 0.000000,},
			{ 0.000000, 0.039655, 0.000000,},
			{ 0.000000, 0.000000, 0.027224,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.071991, 0.143982,},
		R_HEAD = 	{ 0.000000, -0.071991, 0.143982,},
		SGL    = 	{ 0.071991, 0.000000, 0.143982,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.191976, 0.191976, 0.251969,},
		mesh_center = 	{ 0.000000, 0.000000, 0.095988,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.127000, 0.196180,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.710000,
		com = 	{ 0.000000, 0.000000, -0.160636,},
		inertia = 
			{{ 0.017049, 0.000000, 0.000000,},
			{ 0.000000, 0.015188, 0.000000,},
			{ 0.000000, 0.000000, 0.005240,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.027830,},
		R_HUM = 	{ 0.000000, -0.041745, -0.186462,},
		R_HLE = 	{ 0.000000, -0.041745, -0.278302,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.139151, 0.111321, 0.306132,},
		mesh_center = 	{ 0.000000, 0.000000, -0.139151,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.278302,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.620000,
		com = 	{ 0.000000, 0.000000, -0.121511,},
		inertia = 
			{{ 0.008709, 0.000000, 0.000000,},
			{ 0.000000, 0.008029, 0.000000,},
			{ 0.000000, 0.000000, 0.001674,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ -0.053131, 0.000000, -0.265657,},
		R_RSP = 	{ 0.053131, 0.000000, -0.265657,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.106263, 0.079697, 0.292222,},
		mesh_center = 	{ 0.000000, 0.000000, -0.132828,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.265657,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.610000,
		com = 	{ 0.000000, 0.000000, -0.067274,},
		inertia = 
			{{ 0.001744, 0.000000, 0.000000,},
			{ 0.000000, 0.001161, 0.000000,},
			{ 0.000000, 0.000000, 0.000712,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.148507,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.129943, 0.055690, 0.185634,},
		mesh_center = 	{ 0.000000, 0.000000, -0.092817,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.127000, 0.196180,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.710000,
		com = 	{ 0.000000, 0.000000, -0.160636,},
		inertia = 
			{{ 0.017049, 0.000000, 0.000000,},
			{ 0.000000, 0.015188, 0.000000,},
			{ 0.000000, 0.000000, 0.005240,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.027830,},
		L_HUM = 	{ 0.000000, 0.041745, -0.139151,},
		L_HLE = 	{ 0.000000, 0.041745, -0.278302,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.139151, 0.111321, 0.306132,},
		mesh_center = 	{ 0.000000, -0.000000, -0.139151,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.278302,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.620000,
		com = 	{ 0.000000, 0.000000, -0.121511,},
		inertia = 
			{{ 0.008709, 0.000000, 0.000000,},
			{ 0.000000, 0.008029, 0.000000,},
			{ 0.000000, 0.000000, 0.001674,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.053131, 0.000000, -0.265657,},
		L_RSP = 	{ 0.053131, 0.000000, -0.265657,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.106263, 0.079697, 0.292222,},
		mesh_center = 	{ 0.000000, 0.000000, -0.132828,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.265657,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.610000,
		com = 	{ 0.000000, 0.000000, -0.067274,},
		inertia = 
			{{ 0.001744, 0.000000, 0.000000,},
			{ 0.000000, 0.001161, 0.000000,},
			{ 0.000000, 0.000000, 0.000712,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.148507,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.129943, 0.055690, 0.185634,},
		mesh_center = 	{ 0.000000, 0.000000, -0.092817,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}