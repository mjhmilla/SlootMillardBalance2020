return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {65.0},
	 subject_height = {1.63},
	 subject_weight = {97.00},
	 subject_gender = {"male"},
	 subject_pelvisWidth = {0.2730},
	 subject_hipCenterWidth = {0.2080},
	 subject_shoulderCenterWidth = {0.2610},
	 subject_heelAnkleXOffset = {0.0850},
	 subject_heelAnkleZOffset = {0.1010},
	 subject_shoulderNeckZOffset = {0.0780},
	 subject_footWidth = {0.1190},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Thigh_Proximal_R", body = "Thigh_R", point = {-0.079056, 0.000000, 0.039528,},},
	{name = "Thigh_Distal_R", body = "Thigh_R", point = {-0.079056, 0.000000, -0.098821,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.085000, 0.041650, -0.101000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.085000, -0.041650, -0.101000,},},
	{name = "ForeFoot_Medial_R", body = "Foot_R", point = {0.192891, 0.053550, -0.101000,},},
	{name = "ForeFoot_Lateral_R", body = "Foot_R", point = {0.192891, -0.053550, -0.101000,},},
	{name = "Thigh_Proximal_L", body = "Thigh_L", point = {-0.079056, 0.000000, 0.039528,},},
	{name = "Thigh_Distal_L", body = "Thigh_L", point = {-0.079056, 0.000000, -0.098821,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.085000, -0.041650, -0.101000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.085000, 0.041650, -0.101000,},},
	{name = "ForeFoot_Medial_L", body = "Foot_L", point = {0.192891, -0.053550, -0.101000,},},
	{name = "ForeFoot_Lateral_L", body = "Foot_L", point = {0.192891, 0.053550, -0.101000,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.035184, 0.008796, -0.087960,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.035184, 0.008796, -0.087960,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.035184, -0.008796, -0.087960,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.035184, -0.008796, -0.087960,},},
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
		mass   = 10.834900,
		com = 	{ 0.000000, 0.000000, 0.052996,},
		inertia = 
			{{ 0.076256, 0.000000, 0.000000,},
			{ 0.000000, 0.061210, 0.000000,},
			{ 0.000000, 0.000000, 0.069470,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.102308, 0.068205, 0.081846,},
		R_IAS = 	{ 0.102308, -0.068205, 0.081846,},
		L_IPS = 	{ -0.095487, 0.040923, 0.095487,},
		R_IPS = 	{ -0.095487, -0.040923, 0.095487,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.245700, 0.273000, 0.204616,},
		mesh_center = 	{ 0.000000, 0.000000, 0.034103,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.104000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 13.735200,
		com = 	{ 0.000000, 0.000000, -0.161868,},
		inertia = 
			{{ 0.232296, 0.000000, 0.000000,},
			{ 0.000000, 0.232296, 0.000000,},
			{ 0.000000, 0.000000, 0.047645,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.079056, -0.130443,},
		R_FLE = 	{ 0.000000, -0.059292, -0.395282,},
		R_FME = 	{ 0.000000, 0.059292, -0.395282,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.158113, 0.158113, 0.474338,},
		mesh_center = 	{ 0.000000, 0.000000, -0.197641,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.395282,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 4.200100,
		com = 	{ 0.000000, 0.000000, -0.181174,},
		inertia = 
			{{ 0.044966, 0.000000, 0.000000,},
			{ 0.000000, 0.043192, 0.000000,},
			{ 0.000000, 0.000000, 0.007426,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.082446, -0.020611,},
		R_TTC = 	{ 0.020611, 0.000000, -0.061834,},
		R_FAL = 	{ 0.000000, -0.061834, -0.412228,},
		R_TAM = 	{ 0.000000, 0.061834, -0.412228,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.123668, 0.123668, 0.494674,},
		mesh_center = 	{ 0.000000, 0.000000, -0.206114,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.412228,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.328900,
		com = 	{ 0.021686, 0.000000, -0.050500,},
		inertia = 
			{{ 0.005125, 0.000000, 0.000000,},
			{ 0.000000, 0.004658, 0.000000,},
			{ 0.000000, 0.000000, 0.001193,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.169151, 0.072493, -0.072493,},
		R_FM2 = 	{ 0.241644, 0.000000, -0.060411,},
		R_FM5 = 	{ 0.169151, -0.072493, -0.072493,},
		R_FCC = 	{ -0.072493, 0.000000, -0.024164,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.314138, 0.119000, 0.075750,},
		mesh_center = 	{ 0.072069, 0.000000, -0.063125,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.104000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 13.735200,
		com = 	{ 0.000000, 0.000000, -0.161868,},
		inertia = 
			{{ 0.232296, 0.000000, 0.000000,},
			{ 0.000000, 0.232296, 0.000000,},
			{ 0.000000, 0.000000, 0.047645,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.079056, -0.130443,},
		L_FLE = 	{ 0.000000, 0.059292, -0.395282,},
		L_FME = 	{ 0.000000, -0.059292, -0.395282,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.158113, 0.158113, 0.474338,},
		mesh_center = 	{ 0.000000, 0.000000, -0.197641,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.395282,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 4.200100,
		com = 	{ 0.000000, 0.000000, -0.181174,},
		inertia = 
			{{ 0.044966, 0.000000, 0.000000,},
			{ 0.000000, 0.043192, 0.000000,},
			{ 0.000000, 0.000000, 0.007426,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.082446, -0.020611,},
		L_TTC = 	{ 0.020611, 0.000000, -0.061834,},
		L_FAL = 	{ 0.000000, 0.061834, -0.412228,},
		L_TAM = 	{ 0.000000, -0.061834, -0.412228,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.123668, 0.123668, 0.494674,},
		mesh_center = 	{ 0.000000, 0.000000, -0.206114,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.412228,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.328900,
		com = 	{ 0.021686, 0.000000, -0.050500,},
		inertia = 
			{{ 0.005125, 0.000000, 0.000000,},
			{ 0.000000, 0.004658, 0.000000,},
			{ 0.000000, 0.000000, 0.001193,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.169151, -0.072493, -0.072493,},
		L_FM2 = 	{ 0.241644, 0.000000, -0.060411,},
		L_FM5 = 	{ 0.169151, 0.072493, -0.072493,},
		L_FCC = 	{ -0.072493, 0.000000, -0.024164,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.314138, 0.119000, 0.075750,},
		mesh_center = 	{ 0.072069, 0.000000, -0.063125,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.136411,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 15.840100,
		com = 	{ 0.000000, 0.000000, 0.110928,},
		inertia = 
			{{ 0.149804, 0.000000, 0.000000,},
			{ 0.000000, 0.094586, 0.000000,},
			{ 0.000000, 0.000000, 0.141228,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.116000, 0.000000, 0.100880,},
		MAI = 	{ -0.116000, 0.000000, 0.100880,},
		LV1 = 	{ -0.116000, 0.000000, 0.000000,},
		LV3 = 	{ -0.116000, 0.000000, -0.050440,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.232000, 0.273000, 0.302641,},
		mesh_center = 	{ 0.000000, 0.000000, 0.050440,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.201760,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 15.481200,
		com = 	{ 0.000000, 0.000000, 0.111836,},
		inertia = 
			{{ 0.202840, 0.000000, 0.000000,},
			{ 0.000000, 0.081446, 0.000000,},
			{ 0.000000, 0.000000, 0.171980,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.092800, 0.000000, 0.181332,},
		SJN = 	{ 0.097440, 0.000000, 0.135999,},
		TV2 = 	{ -0.092800, 0.000000, 0.135999,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.232000, 0.274050, 0.237998,},
		mesh_center = 	{ 0.000000, 0.000000, 0.090666,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.226665,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 6.731800,
		com = 	{ 0.000000, 0.000000, 0.113661,},
		inertia = 
			{{ 0.031963, 0.000000, 0.000000,},
			{ 0.000000, 0.034545, 0.000000,},
			{ 0.000000, 0.000000, 0.023716,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.068224, 0.136448,},
		R_HEAD = 	{ 0.000000, -0.068224, 0.136448,},
		SGL    = 	{ 0.068224, 0.000000, 0.136448,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.181931, 0.181931, 0.238784,},
		mesh_center = 	{ 0.000000, 0.000000, 0.090965,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.130500, 0.148665,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.628700,
		com = 	{ 0.000000, 0.000000, -0.152231,},
		inertia = 
			{{ 0.014852, 0.000000, 0.000000,},
			{ 0.000000, 0.013231, 0.000000,},
			{ 0.000000, 0.000000, 0.004565,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.026374,},
		R_HUM = 	{ 0.000000, -0.039561, -0.176706,},
		R_HLE = 	{ 0.000000, -0.039561, -0.263740,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.131870, 0.105496, 0.290114,},
		mesh_center = 	{ 0.000000, 0.000000, -0.131870,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.263740,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.571400,
		com = 	{ 0.000000, 0.000000, -0.115153,},
		inertia = 
			{{ 0.007587, 0.000000, 0.000000,},
			{ 0.000000, 0.006994, 0.000000,},
			{ 0.000000, 0.000000, 0.001458,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ -0.050351, 0.000000, -0.251756,},
		R_RSP = 	{ 0.050351, 0.000000, -0.251756,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.100702, 0.075527, 0.276931,},
		mesh_center = 	{ 0.000000, 0.000000, -0.125878,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.251756,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.591700,
		com = 	{ 0.000000, 0.000000, -0.063753,},
		inertia = 
			{{ 0.001519, 0.000000, 0.000000,},
			{ 0.000000, 0.001011, 0.000000,},
			{ 0.000000, 0.000000, 0.000620,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.140736,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.123144, 0.052776, 0.175920,},
		mesh_center = 	{ 0.000000, 0.000000, -0.087960,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.130500, 0.148665,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.628700,
		com = 	{ 0.000000, 0.000000, -0.152231,},
		inertia = 
			{{ 0.014852, 0.000000, 0.000000,},
			{ 0.000000, 0.013231, 0.000000,},
			{ 0.000000, 0.000000, 0.004565,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.026374,},
		L_HUM = 	{ 0.000000, 0.039561, -0.131870,},
		L_HLE = 	{ 0.000000, 0.039561, -0.263740,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.131870, 0.105496, 0.290114,},
		mesh_center = 	{ 0.000000, -0.000000, -0.131870,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.263740,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.571400,
		com = 	{ 0.000000, 0.000000, -0.115153,},
		inertia = 
			{{ 0.007587, 0.000000, 0.000000,},
			{ 0.000000, 0.006994, 0.000000,},
			{ 0.000000, 0.000000, 0.001458,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.050351, 0.000000, -0.251756,},
		L_RSP = 	{ 0.050351, 0.000000, -0.251756,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.100702, 0.075527, 0.276931,},
		mesh_center = 	{ 0.000000, 0.000000, -0.125878,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.251756,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.591700,
		com = 	{ 0.000000, 0.000000, -0.063753,},
		inertia = 
			{{ 0.001519, 0.000000, 0.000000,},
			{ 0.000000, 0.001011, 0.000000,},
			{ 0.000000, 0.000000, 0.000620,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.140736,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.123144, 0.052776, 0.175920,},
		mesh_center = 	{ 0.000000, 0.000000, -0.087960,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}