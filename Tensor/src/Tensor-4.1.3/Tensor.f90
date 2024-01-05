!
!                   _ooOoo_
!                  o8888888o
!                  88" . "88
!                  (| -_- |)
!                  O\  =  /O
!               ____/`---'\____
!             .'  \\|     |//  `.
!            /  \\|||  :  |||//  \
!           /  _||||| -:- |||||-  \
!           |   | \\\  -  /// |   |
!           | \_|  ''\---/''  |   |
!           \  .-\__  `-`  ___/-. /
!         ___`. .'  /--.--\  `. . __
!      ."" '<  `.___\_<|>_/___.'  >'"".
!     | | :  `- \`.;`\ _ /`;.`/ - ` : | |
!     \  \ `-.   \_ __\ /__ _/   .-` /  /
!======`-.____`-.___\_____/___.-`____.-'======
!                   `=---='
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       Buddha blessed , no BUG 
! Report bugs of the package to sj.dong@outlook.com
!

module Tensor_Tools
	use dimension_tools
	use Pointer_Tools
	use Basic_Tools 
	use QuantumNumber_Tools
	use data_base_Tools
	use Tools
	use mpi
	use U1_Tools
	use Parity_Tools
	use memory_type
	implicit none
	private
	type(memory),private::WorkingMemory

	logical::mointer_order_flag=.true.
	character(len=1)::array_character_divider='|'
	character(len=50)::Tensor_Symmetry_type='null'
	character(len=len_of_name)::SVD_U_leg='SVD.Uleg'
	character(len=len_of_name)::SVD_V_leg='SVD.Vleg'
	character(len=len_of_name)::SVD_S_leg1='SVD.s1'
	character(len=len_of_name)::SVD_S_leg2='SVD.s2'
	character(len=len_of_name)::QR_Q_leg='QR.Qleg'
	character(len=len_of_name)::QR_R_leg='QR.Rleg'
	character(len=len_of_name)::LQ_L_leg='LQ.Lleg'
	character(len=len_of_name)::LQ_Q_leg='LQ.Qleg'



	public::Tensor
	type,extends (Dimension) :: Tensor
		type(DataArray),public::Data
	contains
		procedure::deallocatableDimension
		procedure::emptyDimension
		procedure::getTotalDataAll,getTotalDataith,getTotalDatavec
		procedure::writeExternlData,writeExternlData2
		generic,public::write=>writeExternlData2
		procedure::readExternalData
		procedure,public::readData
		generic,public::getTotalData=>getTotalDataAll,getTotalDataith,getTotalDatavec
		procedure::allocateMomeryTensoriType1,allocateMomeryTensoriType2,allocateMomeryTensoriQN
		procedure::allocateMomeryTensoriType3,allocateMomeryTensoriType4,allocateMomeryTensoriQN2
		generic,public::allocateMomery=>allocateMomeryTensoriType1,allocateMomeryTensoriType2,allocateMomeryTensoriQN&
										,allocateMomeryTensoriType3,allocateMomeryTensoriType4,allocateMomeryTensoriQN2
		procedure::allocateMomeryClassTypeiType1,allocateMomeryClassTypeiType2,allocateMomeryClassTypeiQN
		generic,public::allocateMomeryClassType=>allocateMomeryClassTypeiType1,allocateMomeryClassTypeiType2,&
											allocateMomeryClassTypeiQN
		procedure::allocateBlockMomery,allocateBlockMomery2,allocateBlockMomery3,allocateBlockMomery4,allocateBlockMomery5
		generic,public::setBlockMomery=>allocateBlockMomery,allocateBlockMomery2,allocateBlockMomery3&
										,allocateBlockMomery4,allocateBlockMomery5
		procedure::allocateAllBlock1,allocateAllBlock2,allocateAllBlock3,allocateAllBlock4,&
					allocateAllBlock5,allocateAllBlock6,allocateAllBlock7
		generic,public::allocateAllBlock=>allocateAllBlock1,allocateAllBlock2,allocateAllBlock3,&
								allocateAllBlock4,allocateAllBlock5,allocateAllBlock6,allocateAllBlock7
		procedure::allocateAllBlockClassType1,allocateAllBlockClassType3,allocateAllBlockClassType5
		generic,public::allocateAllBlockClassType=>allocateAllBlockClassType1,allocateAllBlockClassType3,&
											allocateAllBlockClassType5
		procedure,public::readOldVersionTensor
		procedure,public::readOldVersionSymTensor
		procedure,public::getType
		procedure,public::getClassType
		procedure,public::getTotalBlock
		procedure,public::getDataFlag
		procedure::getFlagith,getFlagVec,getFlagDimVec,getTensorFlag
		generic,public::getFlag=>getFlagith,getFlagVec,getFlagDimVec,getTensorFlag
		procedure,public::print=>print_as_matrix
		procedure,public::printBlock=>print_as_Block
		procedure::print_info1,print_info2
		generic,public::info=>print_info1,print_info2
		procedure,public::outAllName
		generic,public::getAllName=>outAllName
		procedure,public::getCharacterLen

		procedure,public::asTensor=>SymTensorToTensor
		procedure,public::pick_value=>TensorToSymTensor

		procedure::getValue0i,getValueSymi,getValueVeci,getValueSymVeci,getValueAlli
		procedure::getValue0s,getValueSyms,getValueVecs,getValueSymVecs,getValueAlls
		procedure::getValue0d,getValueSymd,getValueVecd,getValueSymVecd,getValueAlld
		procedure::getValue0c,getValueSymc,getValueVecc,getValueSymVecc,getValueAllc
		procedure::getValue0z,getValueSymz,getValueVecz,getValueSymVecz,getValueAllz
		procedure::getValue0l,getValueSyml,getValueVecl,getValueSymVecl,getValueAlll
		procedure::getValue0a,getValueSyma,getValueVeca,getValueSymVeca,getValueAlla
		procedure::getValue0T,getValueSymT,getValueVecT,getValueSymVecT
		procedure::getithTensor,getvecTensor,getValueAllTensor
		procedure::SymgetithTensor,SymgetvecTensor
		procedure::fermi_element_vec,fermi_element_vec2
		generic,public::ii=>getValue0i,getValueSymi,getValueVeci,getValueSymVeci,getValueAlli
		generic,public::si=>getValue0s,getValueSyms,getValueVecs,getValueSymVecs,getValueAlls
		generic,public::di=>getValue0d,getValueSymd,getValueVecd,getValueSymVecd,getValueAlld
		generic,public::ci=>getValue0c,getValueSymc,getValueVecc,getValueSymVecc,getValueAllc
		generic,public::zi=>getValue0z,getValueSymz,getValueVecz,getValueSymVecz,getValueAllz
		generic,public::li=>getValue0l,getValueSyml,getValueVecl,getValueSymVecl,getValueAlll
		generic,public::ai=>getValue0a,getValueSyma,getValueVeca,getValueSymVeca,getValueAlla
		generic,public::Ti=>getValue0T,getValueSymT,getValueVecT,getValueSymVecT
		generic,public::i=>getithTensor,getvecTensor,getValueAllTensor
		generic,public::Blocki=>SymgetithTensor,SymgetvecTensor
		generic,public::fi=>fermi_element_vec,fermi_element_vec2

		procedure::iwhichindex,swhichindex,dwhichindex,cwhichindex,&
					zwhichindex,awhichindex,awhichindex2
		generic,public::which=>iwhichindex,swhichindex,dwhichindex,cwhichindex,&
								zwhichindex,awhichindex,awhichindex2


		procedure::allocateDimension
		procedure::allocateTensoriType,allocateTensorCha
		procedure::allocateTensoriType2,allocateTensorCha2
		procedure::allocateTensorChaQN,allocateTensoriQN
		procedure::allocateTensoriT,allocateTensoriT2
		generic,public::allocate=>allocateTensoriType,allocateTensorCha,&
								allocateTensoriType2,allocateTensorCha2,&
								allocateTensorChaQN,allocateTensoriQN,&
								allocateTensoriT,allocateTensoriT2
		procedure::allocateTensorClassTpye1,allocateTensorClassTpye2
		procedure::allocateTensorClassTpyeQN,allocateTensorClassTpyeT
		generic,public::allocateClassType=>allocateTensorClassTpye1,allocateTensorClassTpye2,&
							allocateTensorClassTpyeQN,allocateTensorClassTpyeT
		procedure::iPointer1,iPointer2,iPointer3,iPointer4
		procedure::iPointerBlock1,iPointerBlock2,iPointerBlock3,iPointerBlock4
 		procedure::sPointer1,sPointer2,sPointer3,sPointer4
 		procedure::sPointerBlock1,sPointerBlock2,sPointerBlock3,sPointerBlock4
 		procedure::dPointer1,dPointer2,dPointer3,dPointer4
 		procedure::dPointerBlock1,dPointerBlock2,dPointerBlock3,dPointerBlock4
 		procedure::cPointer1,cPointer2,cPointer3,cPointer4
 		procedure::cPointerBlock1,cPointerBlock2,cPointerBlock3,cPointerBlock4
 		procedure::zPointer1,zPointer2,zPointer3,zPointer4
 		procedure::zPointerBlock1,zPointerBlock2,zPointerBlock3,zPointerBlock4
 		procedure::lPointer1,lPointer2,lPointer3,lPointer4
 		procedure::lPointerBlock1,lPointerBlock2,lPointerBlock3,lPointerBlock4
 		procedure::aPointer1,aPointer2,aPointer3,aPointer4
 		procedure::aPointerBlock1,aPointerBlock2,aPointerBlock3,aPointerBlock4
 		procedure::iPointer_vec,sPointer_vec,dPointer_vec,cPointer_vec
 		procedure::zPointer_vec,lPointer_vec,aPointer_vec
 		procedure::iPointerDimBlock2,sPointerDimBlock2,dPointerDimBlock2,cPointerDimBlock2
 		procedure::zPointerDimBlock2,lPointerDimBlock2,aPointerDimBlock2
 		procedure::iPointerDimBlock3,sPointerDimBlock3,dPointerDimBlock3,cPointerDimBlock3
 		procedure::zPointerDimBlock3,lPointerDimBlock3,aPointerDimBlock3



		generic,public::pointer=>iPointer1,iPointer2,iPointer3,iPointer4,&
								iPointerBlock1,iPointerBlock2,iPointerBlock3,iPointerBlock4,&
 								sPointer1,sPointer2,sPointer3,sPointer4,&
 								sPointerBlock1,sPointerBlock2,sPointerBlock3,sPointerBlock4,&
 								dPointer1,dPointer2,dPointer3,dPointer4,&
 								dPointerBlock1,dPointerBlock2,dPointerBlock3,dPointerBlock4,&
 								cPointer1,cPointer2,cPointer3,cPointer4,&
 								cPointerBlock1,cPointerBlock2,cPointerBlock3,cPointerBlock4,&
 								zPointer1,zPointer2,zPointer3,zPointer4,&
 								zPointerBlock1,zPointerBlock2,zPointerBlock3,zPointerBlock4,&
 								lPointer1,lPointer2,lPointer3,lPointer4,&
 								lPointerBlock1,lPointerBlock2,lPointerBlock3,lPointerBlock4,&
 								aPointer1,aPointer2,aPointer3,aPointer4,&
 								aPointerBlock1,aPointerBlock2,aPointerBlock3,aPointerBlock4,&
 								iPointer_vec,sPointer_vec,dPointer_vec,cPointer_vec,&
 								zPointer_vec,lPointer_vec,aPointer_vec,&
 								iPointerDimBlock2,sPointerDimBlock2,dPointerDimBlock2,&
 								cPointerDimBlock2,zPointerDimBlock2,lPointerDimBlock2,&
 								aPointerDimBlock2,&
 								iPointerDimBlock3,sPointerDimBlock3,dPointerDimBlock3,&
 								cPointerDimBlock3,zPointerDimBlock3,lPointerDimBlock3,&
 								aPointerDimBlock3


 		procedure::Tensor1DclassPointer,Tensor2DclassPointer,Tensor3DclassPointer,Tensor4DclassPointer
 		procedure::TensorBlock1DclassPointer,TensorBlock2DclassPointer
 		procedure::TensorBlock3DclassPointer,TensorBlock4DclassPointer
 		procedure::TensorBlockclassPointer
 		procedure::TensorDimBlock2DclassPointer,TensorDimBlock3DclassPointer
 		generic,public::ClassPointer=>Tensor1DclassPointer,Tensor2DclassPointer,Tensor3DclassPointer,&
 									Tensor4DclassPointer,&
 									TensorBlock1DclassPointer,TensorBlock2DclassPointer,&
 									TensorBlock3DclassPointer,TensorBlock4DclassPointer,&
 									TensorBlockclassPointer,&
 									TensorDimBlock2DclassPointer,TensorDimBlock3DclassPointer


 		procedure,public::random=>randomTensorElement


 		procedure::GetBlockDimension
 		generic,public::GetBlockDim=>GetBlockDimension

 		procedure::permutation_vec_routine,permutation_cha_routine,permutation_cha_routine2,permutation_vec_routine2
 		generic,public::permute=>permutation_vec_routine,permutation_cha_routine,permutation_cha_routine2,permutation_vec_routine2
 		procedure::PermuteForWard_ith_routine,PermuteForWard_cha_routine,&
 					PermuteForWard_vec_routine,PermuteForWard_veccha_routine,PermuteForWard_vec_routine2,&
 					PermuteForWard_ith_routine2,PermuteForWard_veccha_routine2,PermuteForWard_cha_routine2
		generic,public::forward=>PermuteForWard_ith_routine,PermuteForWard_cha_routine,&
 					PermuteForWard_vec_routine,PermuteForWard_vec_routine2,PermuteForWard_veccha_routine,&
 					PermuteForWard_ith_routine2,PermuteForWard_veccha_routine2,PermuteForWard_cha_routine2
 		procedure::PermuteBackWard_ith_routine,PermuteBackWard_cha_routine,&
 					PermuteBackWard_vec_routine,PermuteBackWard_veccha_routine,PermuteBackWard_vec_routine2,&
 					PermuteBackWard_ith_routine2,PermuteBackWard_veccha_routine2,PermuteBackWard_cha_routine2
 		generic,public::backward=>PermuteBackWard_ith_routine,PermuteBackWard_cha_routine,&
 					PermuteBackWard_vec_routine,PermuteBackWard_vec_routine2,PermuteBackWard_veccha_routine,&
 					PermuteBackWard_ith_routine2,PermuteBackWard_veccha_routine2,PermuteBackWard_cha_routine2

 		procedure::NotFermiPermuteBackWard_ith_routine,NotFermiPermuteBackWard_vec_routine
 		generic,public::NotFermibackward=>NotFermiPermuteBackWard_ith_routine,NotFermiPermuteBackWard_vec_routine
 		procedure::NotFermiPermuteForWard_ith_routine,NotFermiPermuteForWard_vec_routine
		generic,public::NotFermiforward=>NotFermiPermuteForWard_ith_routine,NotFermiPermuteForWard_vec_routine
		procedure::NotFermipermutation_vec_routine,NotFermipermutation_cha_routine
		procedure::Notfermipermutation_vec_routine2,Notfermipermutation_cha_routine2
		generic,public::NotFermiPermute=>NotFermipermutation_vec_routine,NotFermipermutation_cha_routine,&
										Notfermipermutation_vec_routine2,Notfermipermutation_cha_routine2

 		procedure,public::SymmetryCheck
 		procedure::setZeroValue
 		generic,public::zero=>setZeroValue

 		procedure::contract_vecT,contract_vecA,contract_vecB,contract_ithT,contract_ithA,contract_ithB
 		procedure::contract_ChaVecT,contract_ChaVecA,contract_ChaVecB
 		procedure::contract_ChaT,contract_ChaA,contract_ChaB,contract_Same_nameA
 		procedure::contract_name_ownlegs_routine,contract_ownlegs_routine
 		generic,public::contract=>contract_vecT,contract_vecA,contract_vecB,contract_ithT,contract_ithA,contract_ithB,&
 								  contract_ChaVecT,contract_ChaVecA,contract_ChaVecB,&
 								  contract_ChaT,contract_ChaA,contract_ChaB,contract_Same_nameA,&
 								  contract_name_ownlegs_routine,contract_ownlegs_routine
 		procedure::resetdim1,resetdim2,resetdim3
 		generic,public::resetdim=>resetdim1,resetdim2,resetdim3


 		procedure::subTensorDimSubroutine
 		generic,public::subTensor=>subTensorDimSubroutine

 		procedure::subTensorDeg1DSubroutine,subTensorDeg2DSubroutine,subTensorDeg3DSubroutine
 		generic,public::subDegTensor=>subTensorDeg1DSubroutine,subTensorDeg2DSubroutine,subTensorDeg3DSubroutine

 		procedure::subTensorMaxDimSubroutine
 		generic,public::CutTensorDim=>subTensorMaxDimSubroutine
 		procedure::subTensorMaxDeg1DSubroutine,subTensorMaxDeg2DSubroutine,subTensorMaxDeg3DSubroutine
 		procedure::subTensorCutDeg
 		generic,public::CutTensorDeg=>subTensorMaxDeg1DSubroutine,subTensorMaxDeg2DSubroutine,subTensorMaxDeg3DSubroutine,&
 									subTensorCutDeg


 		procedure::FuseLegs0,FuseLegs1,FuseLegs2
 		generic,public::Fuse=>FuseLegs0,FuseLegs1,FuseLegs2
 		procedure::FuseLegs1Subroutine,FuseLegs2Subroutine,FuseLegs3Subroutine
 		procedure::FuseLegs4Subroutine,FuseLegs5Subroutine,FuseLegs6Subroutine
 		generic,public::FuseTensor=>FuseLegs1Subroutine,FuseLegs2Subroutine,FuseLegs3Subroutine,&
 									FuseLegs4Subroutine,FuseLegs5Subroutine,FuseLegs6Subroutine
 		procedure::SplitTensor1,SplitTensor2,SplitTensor3
 		generic,public::Split=>SplitTensor1,SplitTensor2,SplitTensor3
 		procedure::SplitTensor1Subroutine,SplitTensor2Subroutine,SplitTensor3Subroutine
 		generic,public::SplitTensor=>SplitTensor1Subroutine,SplitTensor2Subroutine,SplitTensor3Subroutine

 		procedure::SVDMatrixNum
 		procedure::SVDNumName,SVDValueName
 		procedure::SVDNumLegs,SVDValueLegs
 		procedure::SVDNumNameMatrixS,SVDValueNameMatrixS
 		procedure::SVDNumLegsMatrixS,SVDValueLegsMatrixS
 		generic,public::SVD=>SVDNumName,SVDValueName,SVDNumNameMatrixS,SVDValueNameMatrixS,&
 							SVDNumLegs,SVDValueLegs,SVDMatrixNum,SVDNumLegsMatrixS,SVDValueLegsMatrixS


 		procedure::eyeQN2,eyeQN3,eyeQN4,eyeQN5,eyeMN1,eyeMN2,eyeMN3,eyeMN4,eye0
 		procedure::eyeQN6,eyeQN7,eyeQN8,eyeQN9
 		generic,public::eye=>eyeQN2,eyeQN3,eyeQN4,eyeQN5,eyeMN1,eyeMN2,eyeMN3,eyeMN4,eye0&
 							,eyeQN6,eyeQN7,eyeQN8,eyeQN9

 		procedure::QRSubroutine
 		generic,public::QRKill=>QRSubroutine

 		procedure::QRSubroutineName,QRSubroutineNumLegs,QRMatrix
 		generic,public::QR=>QRSubroutineName,QRSubroutineNumLegs,QRMatrix

 		procedure::LQSubroutine
 		generic,public::LQkill=>LQSubroutine

 		procedure::LQSubroutineName,LQSubroutineNumLegs,LQmatrix
 		generic,public::LQ=>LQSubroutineName,LQSubroutineNumLegs,LQmatrix

 		procedure::HtransposeTensor,HtransposeTensor2,conjugateTensor
 		generic,public::conjg=>conjugateTensor
 		generic,public::Ntranspose=>HtransposeTensor2
 		generic,public::transpose=>HtransposeTensor

 		procedure::TplusNumSubroutinei,TplusNumSubroutines,TplusNumSubroutined
 		procedure::TplusNumSubroutinec,TplusNumSubroutinez,plusSubroutine
 		generic,public::plus=>TplusNumSubroutinei,TplusNumSubroutines,TplusNumSubroutined,&
 								TplusNumSubroutinec,TplusNumSubroutinez,plusSubroutine
 		procedure::TminusNumSubroutinei,TminusNumSubroutines,TminusNumSubroutined
 		procedure::TminusNumSubroutinec,TminusNumSubroutinez
 		generic,public::minus=>TminusNumSubroutinei,TminusNumSubroutines,TminusNumSubroutined,&
 								TminusNumSubroutinec,TminusNumSubroutinez
 		procedure::NumminusTSubroutinei,NumminusTSubroutines,NumminusTSubroutined
 		procedure::NumminusTSubroutinec,NumminusTSubroutinez				
 		generic,public::minusT=>NumminusTSubroutinei,NumminusTSubroutines,NumminusTSubroutined,&
 								NumminusTSubroutinec,NumminusTSubroutinez					
 		procedure::multiplyi,multiplys,multiplyd,multiplyc,multiplyz
 		generic,public::multiply=>multiplyi,multiplys,multiplyd,multiplyc,multiplyz
 		procedure::dividei,divides,divided,dividec,dividez,TdivideTsubroutine
 		generic,public::divide=>dividei,divides,divided,dividec,dividez,TdivideTsubroutine


 		procedure,public::imax=>maxi
 		procedure,public::smax=>maxs
 		procedure,public::dmax=>maxd
 		procedure,public::cmax=>maxc
 		procedure,public::zmax=>maxz
 		procedure,public::max=>maxT

 		procedure,public::imin=>mini
 		procedure,public::smin=>mins
 		procedure,public::dmin=>mind
 		procedure,public::cmin=>minc
 		procedure,public::zmin=>minz
 		procedure,public::min=>minT

 		procedure,public::isum=>sumi
 		procedure,public::ssum=>sums
 		procedure,public::dsum=>sumd
 		procedure,public::csum=>sumc
 		procedure,public::zsum=>sumz
 		procedure,public::sum=>sumT

 		procedure,public::itrace=>tracei
 		procedure,public::strace=>traces
 		procedure,public::dtrace=>traced
 		procedure,public::ctrace=>tracec
 		procedure,public::ztrace=>tracez
 		procedure,public::trace=>traceT

 		procedure,public::inorm=>normi
 		procedure,public::snorm=>norms
 		procedure,public::dnorm=>normd
 		procedure,public::cnorm=>normc
 		procedure,public::znorm=>normz
 		procedure,public::norm=>normT

 		procedure,public::inorm2=>norm2i
 		procedure,public::snorm2=>norm2s
 		procedure,public::dnorm2=>norm2d
 		procedure,public::cnorm2=>norm2c
 		procedure,public::znorm2=>norm2z
 		procedure,public::norm2=>norm2T

		procedure::set_All_value,set_ith_value0,set_ith_value1,set_vec_value0,set_vec_value1
		procedure::set_some_value0,set_some_value1,set_All_value2
		generic,public::setValue=>set_All_value,set_ith_value0,set_ith_value1,set_vec_value0,&
							set_vec_value1,set_some_value0,set_some_value1,set_All_value2

		procedure::set_tensor_Some_value0,set_tensor_Some_value1,set_tensor_All_value0	
		procedure::set_tensor_value0,set_tensor_vec_value0					
		generic,public::set_Value=>set_tensor_Some_value0,set_tensor_Some_value1,set_tensor_All_value0,&
								set_tensor_value0,set_tensor_vec_value0,set_All_value2	
 								
 		procedure::set_tensor_block_value1,set_tensor_block_value0
 		procedure::set_tensor_blockvec_value0,set_tensor_blockvec_value1
 		generic,public::set_BlockValue=>set_tensor_block_value1,set_tensor_block_value0,&
 										set_tensor_blockvec_value0,set_tensor_blockvec_value1,&
 										set_All_value2

 		procedure::set_tensor_Some_Data,set_tensor_All_data,set_tensor_data,set_tensor_vec_data
 		procedure::set_tensor_block_data,set_tensor_blockvec_data,set_tensor_blockvec_data2
 		procedure::set_tensor_blockvec_data3,set_tensor_All_data2
 		generic,public::setData=>set_tensor_Some_Data,set_tensor_All_data,set_tensor_data,&
 						set_tensor_vec_data,set_tensor_block_data,set_tensor_blockvec_data,&
 						set_tensor_blockvec_data2,set_tensor_blockvec_data3,&
 						set_tensor_All_data2
 		procedure::get_a_Data,get_a_Data_Vec,get_a_Block_Data,get_all_Data
 		procedure::get_Blocki_data,get_Blockvec_data
 		generic,public::getData=>get_a_Data,get_a_Data_Vec,get_a_Block_Data,get_all_Data,&
 								get_Blocki_data,get_Blockvec_data

 		procedure,public::iszero

 		procedure::SubTensor1Func,SubTensor2Func
 		generic,public::getSubTensor=>SubTensor1Func,SubTensor2Func
 		procedure::reOrderToDiag
 		procedure::reOrderTensor
 		procedure::Reverse_Fermi_Rule_specify5,Reverse_Fermi_Rule_specify3
 		generic,public::ReverseFermiArrow=>Reverse_Fermi_Rule_specify5,Reverse_Fermi_Rule_specify3
 		procedure::hermitian_conjugate_Tensor_dimension1
 		procedure::hermitian_conjugate_Tensor_dimension2
 		procedure::hermitian_conjugate_Tensor_dimension3
 		generic,public::hermitian_conjugate_Tensor_dimension=>hermitian_conjugate_Tensor_dimension1,&
 									hermitian_conjugate_Tensor_dimension2,&
 									hermitian_conjugate_Tensor_dimension3
 		procedure,public::setDynamic
		procedure::setType1,setType2
		generic,public::setType=>setType1,setType2
		generic,public::setClassType=>setType1,setType2
		procedure,public::isnan=>isnan0
		procedure,public::linequ
		procedure,public::Solvelinequ
		procedure,public::LLS
		procedure,public::SolveLLS
		procedure::inverseTensor
		generic,public::invTensor=>inverseTensor

		procedure,public::eigTensor
		procedure,public::eig =>eigTensor
		procedure,public::eigRoutine=>eigvalue

		procedure::sortTensor1,sortTensor2,sortTensor3,sortTensor4,sortTensor5
		generic,public::sort=>sortTensor1,sortTensor2,sortTensor3,sortTensor4,sortTensor5

		procedure::enlargeTensor
		generic,public::enlarge=>enlargeTensor

		procedure::store_All_Data,store_a_Tensor
		generic,public::store=>store_All_Data,store_a_Tensor
		procedure::distribute_All_Data,distribute_a_tensor
		generic,public::distribute=>distribute_All_Data,distribute_a_tensor
	end type Tensor


	type(Tensor),target::TMPTensor1,TMPTensor2,TMPproductA,TMPproductB,TMPdirectProductA,TMPdirectProductB,TMPContract1,TMPContract2
	type(Tensor),target::TMPSVD1,TMPSVD2,TMPSVD3,TMPSVD4,eye_A,TMPQR1,TMPQR2,TMPQR3
	type(Tensor),target::SVD_U,SVD_S,SVD_V,TMPpermutation_index
	type(Tensor),target::TMPenlargeTensor
	type(dimension),target::TMPFuseDimension,TMPSplitDimension
	type(Tensor),target::TMPFuseTensor,TMPSplitTensor
	integer,target,allocatable::TMPi1(:),TMPi2(:),SVDCUTDeg(:),SVDPermute(:)
	character(len=Len_of_Name),allocatable::TMPa1(:)
	integer,target,allocatable::product_working_BlockM(:)
	integer,target,allocatable::product_working_BlockN(:)
	integer,target,allocatable::product_working_ResBlockNum(:)
	integer,target,allocatable::reorder_in_QR(:),reorder_in_LQ(:),reorder_in_SVD(:)


	interface Tensor
		procedure constructor0i
		procedure constructor0s
		procedure constructor0d
		procedure constructor0c
		procedure constructor0z
		procedure constructor0l
		procedure constructor_char_scal

		procedure constructor1i
		procedure constructor1s
		procedure constructor1d
		procedure constructor1c
		procedure constructor1z
		procedure constructor1l
		procedure constructor1a

		procedure constructor2i
		procedure constructor2s
		procedure constructor2d
		procedure constructor2c
		procedure constructor2z
		procedure constructor2l
		procedure constructor2a

		procedure constructor3i
		procedure constructor3s
		procedure constructor3d
		procedure constructor3c
		procedure constructor3z
		procedure constructor3l
		procedure constructor3a

		procedure constructor4i
		procedure constructor4s
		procedure constructor4d
		procedure constructor4c
		procedure constructor4z
		procedure constructor4l
		procedure constructor4a
	end interface

	interface
		subroutine tensor_transpose_class_interface(idata,odata,dim_i,plan,rank)
				class(*), intent(in) :: idata(:)
				class(*), intent(out) :: odata(:)
				integer, intent(in) :: dim_i(:), plan(:), rank
		end subroutine tensor_transpose_class_interface
	end interface
	procedure(tensor_transpose_class_interface),public,pointer::tensor_transpose_class=>defaulttensor_transpose_class

	INTERFACE
	  SUBROUTINE checkSymmetryRuleExternal(D1,i1,D2,i2)
			import :: dimension
			type(Dimension),intent(in)::D1,D2
			integer,intent(in)::i1(:),i2(:)
	  END SUBROUTINE checkSymmetryRuleExternal
	END INTERFACE
	procedure(checkSymmetryRuleExternal),public,pointer::checkSymmetryRule=>defaultcheckSymmetryRule





	INTERFACE
	  SUBROUTINE hermitian_conjugate_dimension_interface(dimen,legi)
	  	import :: dimension
		Type(dimension),intent(inout)::dimen
		integer,optional,intent(in)::legi
	  END SUBROUTINE hermitian_conjugate_dimension_interface
	END INTERFACE
	procedure(hermitian_conjugate_dimension_interface),public,pointer::hermitian_conjugate_dimension=>&
															default_hermitian_conjugate_dimension




	interface
		function RuleFuncExternal(dimen,indices)result(Res)
	  		import :: Dimension
	  		logical::Res
			type(Dimension),intent(in)::dimen
			integer,intent(in)::indices(:)
		end function RuleFuncExternal
	end interface
	procedure(RuleFuncExternal),public,pointer::if_Symmetry_Rule=>defaultRule



	interface
		subroutine ifParityExternal(Res,dimen,vec,legi)
			import :: Dimension
			type(Dimension),intent(in)::dimen
			logical,intent(inout)::Res
			integer,intent(in)::vec(:),legi(:)
		end subroutine ifParityExternal
	end interface
	procedure(ifParityExternal),public,pointer::ifParity=>defaultifParity

	interface
		 subroutine QaunNumParityExternal(Res,dimen,ith,jth)
			import :: Dimension
			integer,intent(inout)::Res
			type(Dimension),intent(in)::dimen
			integer,intent(in)::ith,jth
		end subroutine QaunNumParityExternal
	end interface
	procedure(QaunNumParityExternal),public,pointer::QaunNumParity=>defaultQaunNumParity

	interface
		 subroutine SymmetryNewQaunNumInterFace(NewQN,QN)
			real*4,intent(inout)::NewQN
			real*4,intent(in)::QN(:)
		end subroutine SymmetryNewQaunNumInterFace
	end interface
	procedure(SymmetryNewQaunNumInterFace),public,pointer::SymmetryNewQaunNum=>defaultSymmetryNewQaunNum





	public::assignment(=)
	interface assignment(=)
		module procedure TensorToTensor
		module procedure TensorToTensorArrray1
		module procedure TensorToTensorArrray2
		module procedure TensorToTensorArrray3
		module procedure i2T0
		module procedure i2T1
		module procedure i2T2
		module procedure i2T3
		module procedure i2T4
		module procedure T2i0
		module procedure T2i1
		module procedure T2i2
		module procedure T2i3
		module procedure T2i4
		module procedure s2T0
		module procedure s2T1
		module procedure s2T2
		module procedure s2T3
		module procedure s2T4
		module procedure T2s0
		module procedure T2s1
		module procedure T2s2
		module procedure T2s3
		module procedure T2s4
		module procedure d2T0
		module procedure d2T1
		module procedure d2T2
		module procedure d2T3
		module procedure d2T4
		module procedure T2d0
		module procedure T2d1
		module procedure T2d2
		module procedure T2d3
		module procedure T2d4
		module procedure c2T0
		module procedure c2T1
		module procedure c2T2
		module procedure c2T3
		module procedure c2T4
		module procedure T2c0
		module procedure T2c1
		module procedure T2c2
		module procedure T2c3
		module procedure T2c4
		module procedure z2T0
		module procedure z2T1
		module procedure z2T2
		module procedure z2T3
		module procedure z2T4
		module procedure T2z0
		module procedure T2z1
		module procedure T2z2
		module procedure T2z3
		module procedure T2z4
		module procedure l2T0
		module procedure l2T1
		module procedure l2T2
		module procedure l2T3
		module procedure l2T4
		module procedure T2l0
		module procedure T2l1
		module procedure T2l2
		module procedure T2l3
		module procedure T2l4
		module procedure a2T0
		module procedure a2T1
		module procedure a2T2
		module procedure a2T3
		module procedure a2T4
		module procedure T2a0
		module procedure T2a1
		module procedure T2a2
		module procedure T2a3
		module procedure T2a4
	end interface



	public::operator(.p.)
	interface operator(.p.)
		module procedure permutation_vec
		module procedure permutation_cha
	end interface

	public::operator(.pf.)
	interface operator(.pf.)
		module procedure PermuteForWard_ith
		module procedure PermuteForWard_cha
		module procedure PermuteForWard_vec
		module procedure PermuteForWard_veccha
	end interface

	public::operator(.pb.)
	interface operator(.pb.)
		module procedure PermuteBackWard_ith
		module procedure PermuteBackWard_cha
		module procedure PermuteBackWard_vec
		module procedure PermuteBackWard_veccha
	end interface

	public::operator(.npf.)
	interface operator(.npf.)
		module procedure NotFermiPermuteForWard_ith
		module procedure NotFermiPermuteForWard_vec
	end interface

	public::operator(.npb.)
	interface operator(.npb.)
		module procedure NotFermiPermuteBackWard_ith
		module procedure NotFermiPermuteBackWard_vec
	end interface

	public::operator(.np.)
	interface operator(.np.)
		module procedure NotFermipermutation_vec
		module procedure NotFermipermutation_cha
	end interface


	public::contract
	interface contract
		module procedure contract_vec
		module procedure contract_ith
		module procedure contract_ChaVec
		module procedure contract_Cha
		module procedure contract_Same_name
		module procedure contract_ownlegs
		module procedure contract_name_ownlegs
	end interface


	interface tensor_transpose
		procedure tensor_transpose_real4
		procedure tensor_transpose_real8
		procedure tensor_transpose_complex8
		procedure tensor_transpose_complex16
		procedure tensor_transpose_int
		procedure tensor_transpose_logical
		procedure tensor_transpose_string
	end interface


	public::Reverse_Fermi_Rule
	interface Reverse_Fermi_Rule
		module procedure Reverse_Fermi_Rule_specify2
		module procedure Reverse_Fermi_Rule_specify3
		module procedure Reverse_Fermi_Rule_specify4
		module procedure Reverse_Fermi_Rule_specify5
		module procedure Reverse_Fermi_Rule_specify6
		module procedure Reverse_Fermi_Rule_specify7
		module procedure Reverse_Fermi_Rule1
	end interface

	public::ReverseFermiArrow
	interface ReverseFermiArrow
		module procedure Reverse_Fermi_Rule_specify2
		module procedure Reverse_Fermi_Rule_specify3
		module procedure Reverse_Fermi_Rule_specify4
		module procedure Reverse_Fermi_Rule_specify5
		module procedure Reverse_Fermi_Rule_specify6
		module procedure Reverse_Fermi_Rule_specify7
		module procedure Reverse_Fermi_Rule1
	end interface

	public::expm
	interface expm
		module procedure  expmMatrix1
	end interface

	public::isnan
	interface isnan
		module procedure  isnan
	end interface

	public::operator(.H.)
	interface operator(.H.)
		module procedure Htranspose
	end interface
	public::operator(.Hn.)
	interface operator(.Hn.)
		module procedure Htranspose2
	end interface
	public::operator(.con.)
	interface operator(.con.)
		module procedure conjugate! conjugate
	end interface

	public::conjg
	interface conjg
		module procedure conjugate! conjugate
	end interface

	public::abs
	interface abs
		module procedure absTensor
	end interface

	public::dabs
	interface dabs
		module procedure dabsTensor
	end interface

	public::eye
	interface eye
		module procedure eyeQNTensor1
		module procedure eyeQNTensor2
		module procedure eyeQNTensor3
		module procedure eyeQNTensor4
		module procedure eyeQNTensor5
		module procedure eyeQNTensor6
		module procedure eyeMNTensor1
		module procedure eyeMNTensor2
		module procedure eyeMNTensor3
	end interface

	public::operator(.x.)
	interface operator(.x.)! dot product conjugating the first vector,The Tensor will be regard as a vector
		module procedure product_cdotT
	end interface
	public::operator(.ix.)
	interface operator(.ix.)
		module procedure product_cdoti
	end interface
	public::operator(.sx.)
	interface operator(.sx.)
		module procedure product_cdots
	end interface
	public::operator(.dx.)
	interface operator(.dx.)
		module procedure product_cdotd
	end interface
	public::operator(.cx.)
	interface operator(.cx.)
		module procedure product_cdotc
	end interface
	public::operator(.zx.)
	interface operator(.zx.)
		module procedure product_cdotz
	end interface

	public::operator(.dot.)
	interface operator(.dot.)! dot product DO NOT conjugating the first vector,The Tensor will be regard as a vector
		module procedure product_dotT
	end interface
	public::operator(.idot.)
	interface operator(.idot.)
		module procedure product_doti
	end interface
	public::operator(.sdot.)
	interface operator(.sdot.)
		module procedure product_dots
	end interface
	public::operator(.ddot.)
	interface operator(.ddot.)
		module procedure product_dotd
	end interface
	public::operator(.cdot.)
	interface operator(.cdot.)
		module procedure product_dotc
	end interface
	public::operator(.zdot.)
	interface operator(.zdot.)
		module procedure product_dotz
	end interface

	public::operator(+)
	interface operator(+)
		module procedure plusTensor
		module procedure TplusNumi
		module procedure TplusNums
		module procedure TplusNumd
		module procedure TplusNumc
		module procedure TplusNumz
		module procedure TplusNuma
		module procedure NumplusTi
		module procedure NumplusTs
		module procedure NumplusTd
		module procedure NumplusTc
		module procedure NumplusTz
		module procedure NumplusTa
	end interface

	public::operator(-)
	interface operator(-)
		module procedure minusTensor
		module procedure TMinusNumi
		module procedure TMinusNums
		module procedure TMinusNumd
		module procedure TMinusNumc
		module procedure TMinusNumz
		module procedure NumMinusTi
		module procedure NumMinusTs
		module procedure NumMinusTd
		module procedure NumMinusTc
		module procedure NumMinusTz
	end interface

	public::operator(*)
	interface operator(*)
		module procedure ProductTensor
		module procedure TmultiplyNumi
		module procedure TmultiplyNums
		module procedure TmultiplyNumd
		module procedure TmultiplyNumc
		module procedure TmultiplyNumz
		module procedure NummultiplyTi
		module procedure NummultiplyTs
		module procedure NummultiplyTd
		module procedure NummultiplyTc
		module procedure NummultiplyTz
	end interface

	public::operator(/)
	interface operator(/)
		module procedure TdivideNumi
		module procedure TdivideNums
		module procedure TdivideNumd
		module procedure TdivideNumc
		module procedure TdivideNumz
		module procedure TdivideT
		module procedure NumdivideTi
		module procedure NumdivideTs
		module procedure NumdivideTd
		module procedure NumdivideTc
		module procedure NumdivideTz
	end interface

	public::operator(.kron.)!direct Product,support any rank tensor and keep their TensorName,see more in help/operator
	interface operator(.kron.)
		module procedure directProductTensor
	end interface


	public::set_symmetry,ProductTensorRoutine,ProductTensorCase12,MatrixDim
	public::identityMatrix,out_Tensor_Symmetry_type
	public::set_array_character_divider,get_array_character_divider
	public::send_Tensor,BCAST_Tensor,ALLREDUCE_Tensor,REDUCE_Tensor
	public::MPI_send_Tensor,MPI_BCAST_Tensor
	interface MPI_send_Tensor
		module procedure send_Tensor
	end interface
	interface MPI_BCAST_Tensor
		module procedure BCAST_Tensor
	end interface

	public::MPI_SUM_Tensor
	interface MPI_SUM_Tensor
		module procedure MPI_SUM_Tensor1
		module procedure MPI_SUM_Tensor2
	end interface

	public::MPI_MAX_Tensor
	interface MPI_MAX_Tensor
		module procedure MPI_MAX_Tensor1
		module procedure MPI_MAX_Tensor2
	end interface

	public::MPI_MIN_Tensor
	interface MPI_MIN_Tensor
		module procedure MPI_MIN_Tensor1
		module procedure MPI_MIN_Tensor2
	end interface

	public::int
	interface int
		module procedure intTensor
	end interface

	public::char
	interface char
		module procedure charTensor
	end interface

	public::real
	interface real
		module procedure realTensor
	end interface

	public::dble
	interface dble
		module procedure dbleTensor
	end interface

	public::dreal
	interface dreal
		module procedure dbleTensor
	end interface

	public::imag
	interface imag
		module procedure imagTensor
	end interface

	public::dimag
	interface dimag
		module procedure dimagTensor
	end interface

	public::cmplx
	interface cmplx
		module procedure cmplxTensor
		module procedure cmplxTensor2
	end interface

	public::dcmplx
	interface dcmplx
		module procedure dcmplxTensor
		module procedure dcmplxTensor2
	end interface

	public::eSquare
	interface eSquare
		module procedure elementSquare
	end interface

	!eSquare2: x+iy -->  x*x + i y*y

	public::eSquare2
	interface eSquare2
		module procedure elementSquare2
	end interface

	public::eNorm2
	interface eNorm2
		module procedure elementnorm2
	end interface

	public::eSqrt
	interface eSqrt
		module procedure elementSqrt
	end interface

	public::eProduct
	interface eProduct
		module procedure elementProduct
	end interface

	public::edivide
	interface edivide
		module procedure elementdivide
	end interface

	public::dsum
	interface dsum
		module procedure sumd
	end interface

	public::isum
	interface isum
		module procedure sumi
	end interface

	public::ssum
	interface ssum
		module procedure sums
	end interface

	public::csum
	interface csum
		module procedure sumc
	end interface

	public::zsum
	interface zsum
		module procedure sumz
	end interface

	

	public::operator(.eq.)
	interface operator(.eq.)
		module procedure equal_of_Tensor
		module procedure T_eq_i
		module procedure T_eq_s
		module procedure T_eq_d
		module procedure i_eq_T
		module procedure s_eq_T
		module procedure d_eq_T
	end interface

	public::operator(.equ.)
	interface operator(.equ.)
		module procedure equ_of_Tensor
		module procedure T_equ_i
		module procedure T_equ_s
		module procedure T_equ_d
		module procedure i_equ_T
		module procedure s_equ_T
		module procedure d_equ_T
	end interface
	
	public::operator(.gt.)
	interface operator(.gt.)
		module procedure gt_of_Tensor
		module procedure T_gt_i
		module procedure T_gt_s
		module procedure T_gt_d
		module procedure i_gt_T
		module procedure s_gt_T
		module procedure d_gt_T
	end interface
	
	public::operator(.ge.)
	interface operator(.ge.)
		module procedure ge_of_Tensor
		module procedure T_ge_i
		module procedure T_ge_s
		module procedure T_ge_d
		module procedure i_ge_T
		module procedure s_ge_T
		module procedure d_ge_T
	end interface
	
	public::operator(.lt.)
	interface operator(.lt.)
		module procedure lt_of_Tensor
		module procedure T_lt_i
		module procedure T_lt_s
		module procedure T_lt_d
		module procedure i_lt_T
		module procedure s_lt_T
		module procedure d_lt_T
	end interface
	
	public::operator(.le.)
	interface operator(.le.)
		module procedure le_of_Tensor
		module procedure T_le_i
		module procedure T_le_s
		module procedure T_le_d
		module procedure i_le_T
		module procedure s_le_T
		module procedure d_le_T
	end interface
	public::operator(.ne.)
	interface operator(.ne.)
		module procedure T_ne_i
		module procedure T_ne_s
		module procedure T_ne_d
		module procedure i_ne_T
		module procedure s_ne_T
		module procedure d_ne_T
	end interface

	public::operator(.nequ.)
	interface operator(.nequ.)
		module procedure T_nequ_i
		module procedure T_nequ_s
		module procedure T_nequ_d
		module procedure i_nequ_T
		module procedure s_nequ_T
		module procedure d_nequ_T
	end interface

	public::operator(.inv.)
	interface operator(.inv.)
		module procedure inverse
		module procedure inverseTen
	end interface

	public::PasteTensor
	interface PasteTensor
		module procedure pasteTensorFunc
	end interface

	public::addTensorcol
	interface addTensorcol
		module procedure combinationColFunc
	end interface

	public::addTensorRow
	interface addTensorRow
		module procedure combinationrowFunc
	end interface
	

	public::writemess
	interface writemess
		module procedure writemess_Tensor
		module procedure writemess_Tensor_form
	end interface
	
	public::pauli_matrix
	public::set_SVD_leg_name_info
	public::set_QR_leg_name_info
	public::set_LQ_leg_name_info
	


contains

	subroutine defaulttensor_transpose_class(idata, odata, dim_i, plan, rank)
		class(*), intent(in) :: idata(:)
		class(*), intent(out) :: odata(:)
		integer, intent(in) :: dim_i(:), plan(:), rank
		integer, allocatable :: ldi(:), ldo(:), dim_o(:), ldi_i(:)
		integer :: i, ii, io, active
		integer, allocatable :: index(:)
		class(*),pointer::op,ip

		allocate(ldi(rank))
		allocate(ldo(rank))
		allocate(dim_o(rank))
		allocate(ldi_i(rank))
		allocate(index(rank))

		do i=1, rank
			dim_o(i) = dim_i(plan(i))
		end do
		ldo(1) = 1
		do i=1, rank- 1
			ldo(i+1) = ldo(i)*dim_o(i)
		end do
		ldi_i(1) = 1
		do i=1, rank - 1
			ldi_i(i+1) = ldi_i(i)*dim_i(i)
		end do
		do i=1, rank
			ldi(i) = ldi_i(plan(i))
		end do

		do i=1, rank
			index(i) = 1
		end do
		ii = 1
		io = 1
		do
			call ClassPointer1DFunc(odata,io,op)
			call ClassPointer1DFunc(idata,ii,ip)
			call CopyData(op,ip)
			!odata(io) = idata(ii)
			active = 1
			index(active) = index(active) + 1
			ii = ii + ldi(active)
			io = io + ldo(active)
			do while (index(active) > dim_o(active))
				index(active) = 1
				ii = ii - dim_o(active) * ldi(active)
				io = io - dim_o(active) * ldo(active)
				if (active == rank) then
					return
				end if
				active = active + 1
				index(active) = index(active) + 1
				ii = ii + ldi(active)
				io = io + ldo(active)
			end do
		end do
	end subroutine

	subroutine set_SVD_leg_name_info(U,S1,S2,V)
		character(len=*),intent(in)::U,S1,S2,V
		call writemess('Set the name info of the SVD output:')
		call writemess('U tensor: '+U)
		call writemess('S tensor: '+S1+', '+S2)
		call writemess('V tensor: '+V)
		SVD_U_leg=U
		SVD_V_leg=V
		SVD_S_leg1=S1
		SVD_S_leg2=S2
		return
	end subroutine
	subroutine set_QR_leg_name_info(Q,R)
		character(len=*),intent(in)::Q,R
		call writemess('Set the name info of the QR output:')
		call writemess('Q tensor: '+Q)
		call writemess('R tensor: '+R)
		QR_Q_leg=Q
		QR_R_leg=R
		return
	end subroutine
	subroutine set_LQ_leg_name_info(L,Q)
		character(len=*),intent(in)::L,Q
		call writemess('Set the name info of the QR output:')
		call writemess('L tensor: '+L)
		call writemess('Q tensor: '+Q)
		LQ_L_leg=L
		LQ_Q_leg=Q
		return
	end subroutine

#include "Tensor_include/eigenvaule.f90"

#include "Tensor_include/linear_equation.f90"

#include "Tensor_include/linear_equation2.f90"

#include "Tensor_include/inverse.f90"

#include "Tensor_include/getData.f90"

#include "Tensor_include/SetData.f90"

#include "Tensor_include/pointer_and_other_func.f90"

#include "Tensor_include/permutation.f90"

#include "Tensor_include/product.f90"

#include "Tensor_include/fermi_arrow.f90"

#include "Tensor_include/contract.f90"

#include "Tensor_include/SetDimension.f90"

#include "Tensor_include/reshape.f90"

#include "Tensor_include/fuse_split.f90"

#include "Tensor_include/SVD.f90"

#include "Tensor_include/eye.f90"

#include "Tensor_include/QR.f90"

#include "Tensor_include/LQ.f90"

#include "Tensor_include/basic_operation.f90"

#include "Tensor_include/assignment.f90"

#include "Tensor_include/getValue.f90"

#include "Tensor_include/constructor.f90"

#include "Tensor_include/OtherFunction.f90"

#include "Tensor_include/subTensor.f90"

#include "Tensor_include/set_value.f90"

#include "Tensor_include/which.f90"

#include "Tensor_include/logicalFunc.f90"

#include "Tensor_include/FuncOnElement.f90"

#include "Tensor_include/sort.f90"

#include "Tensor_include/paste.f90"

#include "Tensor_include/enlarge.f90"

	function get_array_character_divider()
		character(len=1)::get_array_character_divider
		get_array_character_divider=array_character_divider
		return
	end function

	subroutine set_array_character_divider(cha)
		character(len=1),intent(in)::cha
		array_character_divider=cha
		return
	end subroutine

	function out_Tensor_Symmetry_type()
		character(len=50)::out_Tensor_Symmetry_type
		out_Tensor_Symmetry_type=Tensor_Symmetry_type
	end function

	!*********************************************************
	!
	! default symmetry subroutine
	!
	!*********************************************************

	subroutine SymmetryCheck(T)
		class(Tensor),intent(in)::T
		integer::i,rank
		logical::checkFlag
		integer,allocatable::Tdim(:),Maxdim(:)
		integer,pointer::dim(:),si(:)
		if(.not.T%getFlag())return
		if(.not.T%getSymmetryFlag())return
		rank=T%GetRank()
		allocate(TDim(rank))
		allocate(Maxdim(rank))
		Maxdim=T%dim()
		checkFlag=.true.
		do i=1,T%getTotalBlock()
			if(T%getFlag(i))then
				call IndesToaddress(Maxdim,Tdim,i)	
				checkFlag=if_symmetry_rule(T%Dimension,TDim)
				if(.not.checkFlag)then
					call writemess(' The Tensor is not a Symmetry Tensor')
					call writemess('Tdim:')
					call writemess(Tdim)
					call error_stop
				end if
			end if
		end do
		call T%pointDim(dim)
		call T%Data%pointStarti(si)
		if(product(dim).ne.size(si))then
			call writemess(' Internal error in the symmetry tensor')
			call error_stop
		end if
		call T%Data%pointEndi(si)
		if(product(dim).ne.size(si))then
			call writemess(' Internal error in the symmetry tensor')
			call error_stop
		end if
		return
	end subroutine


	subroutine error_mess()
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call writemess('@@@                 ERROR                   @@@  ')
		call writemess('@@@ You have not set the symmetry group yet @@@  ')
		call writemess('@@@                                         @@@  ')
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call writemess('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  ')
		call error_stop
	end subroutine


	subroutine defaultcheckSymmetryRule(D1,i1,D2,i2)
		type(Dimension),intent(in)::D1,D2
		integer,intent(in)::i1(:),i2(:)
	end subroutine


	subroutine default_hermitian_conjugate_dimension(dimen,legi)
		Type(dimension),intent(inout)::dimen
			integer,optional,intent(in)::legi
		call error_mess()
	end subroutine


	function defaultRule(dimen,indices)result(Res)
  		logical::Res
		type(Dimension),intent(in)::dimen
		integer,intent(in)::indices(:)
		call error_mess()
	end function

	subroutine defaultifParity(Res,dimen,vec,legi)
		logical,intent(inout)::Res
		type(Dimension),intent(in)::dimen
		integer,intent(in)::vec(:),legi(:)
		call error_mess()
	end subroutine

	subroutine defaultQaunNumParity(REs,dimen,ith,jth)
		integer,intent(inout)::Res
		type(Dimension),intent(in)::dimen
		integer,intent(in)::ith,jth
		call error_mess()
	end subroutine

	subroutine defaultSymmetryNewQaunNum(NewQN,QN)
		real*4,intent(inout)::NewQN
		real*4,intent(in)::QN(:)
		call error_mess()
	end subroutine 

	subroutine set_symmetry(group)
		character(len=*),intent(in)::group
		call writemess('   ')
		call writemess(' Set the symmetry group as '+(' '+group))
		call writemess('   ')
		if(group.equ.'U(1)')then
			checkSymmetryRule=>checkU1Rule
			if_Symmetry_Rule=>U1RuleFunc
			ifParity=>U1ifParity
			QaunNumParity=>U1Parity
			SymmetryNewQaunNum=>U1NewQN
			hermitian_conjugate_dimension=>hermitian_conjugate_U1_dimension
			call set_SetDimensionRule_function(hermitian_conjugate_U1_dimension)
			Tensor_Symmetry_type='U(1)'
		else if(group.equ.'Parity')then
			checkSymmetryRule=>checkParityRule
			if_Symmetry_Rule=>ParityRule
			ifParity=>ParityifParity
			QaunNumParity=>ParityParity
			SymmetryNewQaunNum=>ParityNewQN
			hermitian_conjugate_dimension=>hermitian_conjugate_parity_dimension
			call set_SetDimensionRule_function(hermitian_conjugate_parity_dimension)
			Tensor_Symmetry_type='Parity'
		else if(group.equ.'null')then
			call writemess(' DO not set the symmetry group')
		else
			call writemess(' Only support the symmetry group of U(1) or Parity')
			call error_stop
		end if
		return
	end subroutine


	!*********************************************************
	!
	! assignment
	!
	!*********************************************************

	subroutine TensorToTensor(A,B)
		type(Tensor),intent(inout)::A
		type(Tensor),intent(in)::B
		if(.not.B%getFlag())then
			call A%empty()
			return
		end if
		A%Dimension=B%dimension
		A%Data=B%Data
		return
	end subroutine

	subroutine TensorToTensorArrray1(A,B)
		type(Tensor),intent(inout)::A(:)
		type(Tensor),intent(in)::B
		integer::i
		if(.not.B%getFlag())then
			do i=1,size(A)
				call A(i)%empty()
			end do
			return
		end if
		do i=1,size(A)
			call TensorToTensor(A(i),B)
		end do
		return
	end subroutine

	subroutine TensorToTensorArrray2(A,B)
		type(Tensor),intent(inout)::A(:)
		type(Tensor),intent(in)::B(:)
		integer::i
		if(size(A).lt.size(B))then
			call writemess('ERROR in assignment of array for tensors',-1)
			call writemess('size(A)<size(B)',-1)
			call writemess('size(A)='+size(A),-1)
			call writemess('size(B)='+size(B),-1)
			call error_Stop
		end if
		do i=1,size(B)
			call TensorToTensor(A(i),B(i))
		end do
		return
	end subroutine
	subroutine TensorToTensorArrray3(A,B)
		type(Tensor),intent(inout)::A
		type(Tensor),intent(in)::B(:)
		integer::i
		if(size(B).ne.1)then
			call writemess('ERROR in assignment of array for tensors',-1)
			call writemess('size(A)<size(B)',-1)
			call writemess('size(A)=1',-1)
			call writemess('size(B)='+size(B),-1)
			call error_Stop
		end if
		call TensorToTensor(A,B(1))
		return
	end subroutine

	!*********************************************************
	!
	!   Setting value
	!
	!*********************************************************

	subroutine SymTensorToTensor(SymT,Res)
		class(Tensor),intent(in)::SymT
		type(Tensor),intent(inout)::Res
		type(Dimension)::Dimen
		integer,pointer::deg(:),inde(:),SymTdim(:),minmaxindex(:,:),ip(:)
		integer::Degi,i,rank,ii,dim1,dim2
		logical::goon
		integer,pointer::iResp(:),iSymTp(:)
		real*4,pointer::sResp(:),sSymTp(:)
		real*8,pointer::dResp(:),dSymTp(:)
		complex*8,pointer::cResp(:),cSymTp(:)
		complex*16,pointer::zResp(:),zSymTp(:)
		logical,pointer::lResp(:),lSymTp(:)
		character(len=SymT%Data%DataCharacterLen),pointer::aResp(:),aSymTp(:)
		integer::iType
		class(*),pointer::mold,p,clsResp(:),clsSymTp(:)
		if(SymT%getTotalBlock().eq.1)then
			iType=SymT%getType()
			call ClassPointer1DFunc(SymT%Data%ClassData,1,mold)
			if((iType.ge.0).and.(iType.le.8))then
				call Res%allocate([1],SymT%getType(),SymT%getCharacterLen())
				call ClassPointer1DFunc(Res%Data%ClassData,1,p)
			else
				call Res%Data%allocateClassType(1,mold,SymT%getType())
				Res%Dimension=[1]
				call ClassPointer1DFunc(Res%Data%ClassData,1,p)
			end if
			call CopyData(p,mold)
			return
		end if
		if(SymT%getTotalBlock().eq.0)then
			iType=SymT%getType()
			if((iType.ge.0).and.(iType.le.8))then
				call Res%allocate([1],SymT%getType(),SymT%getCharacterLen())
			else
				call ClassPointer1DFunc(SymT%Data%ClassData,1,mold)
				call Res%Data%allocateClassType(1,mold,SymT%getType())
				Res%Dimension=[1]
			end if
			call Res%zero()
			return
		end if

		rank=SymT%getRank()
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+rank)
		call WorkingMemory%get_memory(deg,rank)

		
		call Res%empty()
		do i=1,rank
			call SymT%pointDeg(ip,i)
			deg(i)=sum(ip)
		end do
		if(product(deg).eq.0)then
			call WorkingMemory%free()
			return
		end if

		call Res%allocate(deg,SymT%getType())
		

		call WorkingMemory%get_memory(inde,rank)
		call WorkingMemory%get_memory(minmaxindex,2,rank)

		call SymT%pointDim(SymTdim)
		call Res%zero()
		select case(SymT%getType())
			case(1)
				call Res%pointer(iResp)
				do i=1,SymT%getTotalBlock()
					if(SymT%Data%getFlag(i))then
						call IndesToaddress(SymTdim,inde,i)
						call SymT%pointer(iSymTp,i)
						minmaxindex=SymT%NonSymIndex(inde)
						call iStoreData(iResp,iSymTp,minmaxindex(1,:),minmaxindex(2,:),deg)
					end if
				end do
			case(2)
				call Res%pointer(sResp)
				do i=1,SymT%getTotalBlock()
					if(SymT%Data%getFlag(i))then
						call IndesToaddress(SymTdim,inde,i)
						call SymT%pointer(sSymTp,i)
						minmaxindex=SymT%NonSymIndex(inde)
						call sStoreData(sResp,sSymTp,minmaxindex(1,:),minmaxindex(2,:),deg)
					end if
				end do
			case(3)
				call Res%pointer(dResp)
				do i=1,SymT%getTotalBlock()
					if(SymT%Data%getFlag(i))then
						call IndesToaddress(SymTdim,inde,i)
						call SymT%pointer(dSymTp,i)
						minmaxindex=SymT%NonSymIndex(inde)
						call dStoreData(dResp,dSymTp,minmaxindex(1,:),minmaxindex(2,:),deg)
					end if
				end do
			case(4)
				call Res%pointer(cResp)
				do i=1,SymT%getTotalBlock()
					if(SymT%Data%getFlag(i))then
						call IndesToaddress(SymTdim,inde,i)
						call SymT%pointer(cSymTp,i)
						minmaxindex=SymT%NonSymIndex(inde)
						call cStoreData(cResp,cSymTp,minmaxindex(1,:),minmaxindex(2,:),deg)
					end if
				end do
			case(5)
				call Res%pointer(zResp)
				do i=1,SymT%getTotalBlock()
					if(SymT%Data%getFlag(i))then
						call IndesToaddress(SymTdim,inde,i)
						call SymT%pointer(zSymTp,i)
						minmaxindex=SymT%NonSymIndex(inde)
						call zStoreData(zResp,zSymTp,minmaxindex(1,:),minmaxindex(2,:),deg)
					end if
				end do
			case(6)
				call Res%pointer(lResp)
				do i=1,SymT%getTotalBlock()
					if(SymT%Data%getFlag(i))then
						call IndesToaddress(SymTdim,inde,i)
						call SymT%pointer(lSymTp,i)
						minmaxindex=SymT%NonSymIndex(inde)
						call lStoreData(lResp,lSymTp,minmaxindex(1,:),minmaxindex(2,:),deg)
					end if
				end do
			case(7)
				call Res%pointer(aResp)
				do i=1,SymT%getTotalBlock()
					if(SymT%Data%getFlag(i))then
						call IndesToaddress(SymTdim,inde,i)
						call SymT%pointer(aSymTp,i)
						minmaxindex=SymT%NonSymIndex(inde)
						call aStoreData(aResp,aSymTp,minmaxindex(1,:),minmaxindex(2,:),deg)
					end if
				end do
			case default
				call Res%ClassPointer(clsResp)
				do i=1,SymT%getTotalBlock()
					if(SymT%Data%getFlag(i))then
						call IndesToaddress(SymTdim,inde,i)
						call SymT%ClassPointer(clsSymTp,i)
						minmaxindex=SymT%NonSymIndex(inde)
						call ClassStoreData(clsResp,clsSymTp,minmaxindex(1,:),minmaxindex(2,:),deg)
					end if
				end do
		end select

			
		if(.not.SymT%getNameFlag())then
			call WorkingMemory%free()
			return
		end if

		do i=1,rank
			call Res%setName(i,SymT%getName(i))
		end do
		call WorkingMemory%free()
		return
	end subroutine

	subroutine TensorToSymTensor(Res,T)
		class(Tensor),intent(inout)::Res
		type(Tensor),intent(in)::T
		integer,allocatable::inde(:),alldim(:),minmaxindex(:,:),indexstart(:),indexend(:)
		integer,pointer::Tdim(:)
		integer::Degi,i,rank,classtype,ii,dim1,dim2
		logical::goon
		integer,pointer::iRp(:),iTp(:)
		real*4,pointer::sRp(:),sTp(:)
		real*8,pointer::dRp(:),dTp(:)
		complex*8,pointer::cRp(:),cTp(:)
		complex*16,pointer::zRp(:),zTp(:)
		logical,pointer::lRp(:),lTp(:)
		character(len=T%Data%DataCharacterLen),pointer::aRp(:),aTp(:)
		class(*),pointer::ClsTp(:),clsRp(:)
		type(Tensor)::temp1,temp2
		if(.not.Res%getFlag())then
			call writemess("ERROR in SymTensor=Tensor,the SymTensor should be allocate ")
			call error_stop()
		endif
		rank=Res%getRank()
		classtype=T%getType()
		if(classtype.ne.Res%getType())then
			call writemess("ERROR in SymTensor=Tensor, data type do not match")
			call writemess('dataype'+classtype+','+Res%getType())
			call error_stop()
		end if

		allocate(inde(rank))
		allocate(alldim(rank))
		allocate(indexstart(rank))
		allocate(indexend(rank))
		allocate(minmaxindex(2,rank))
		alldim=Res%dim()
		indexstart=1
		call T%pointDim(TDim)
		select case(T%getType())
			case (1)
				call T%pointer(iTp)
				do i=1,Res%getTotalBlock()
					if(Res%getFlag(i))then
						call IndesToaddress(alldim,inde,i)
						minmaxindex=Res%NonSymIndex(inde)
						indexend=Res%getBlockdim(inde)
						!!inoutT(indexstart:indexend)=inoutT(minindeT:maxTinde)
						call Res%pointer(iRp,i)
						call iStoreSomeData(iRp,indexend,iTp,TDim,indexstart,indexend,minmaxindex(1,:),minmaxindex(2,:))
					end if
				end do
			case (2)
				call T%pointer(sTp)
				do i=1,Res%getTotalBlock()
					if(Res%getFlag(i))then
						call IndesToaddress(alldim,inde,i)
						minmaxindex=Res%NonSymIndex(inde)
						indexend=Res%getBlockdim(inde)
						!!inoutT(indexstart:indexend)=inoutT(minindeT:maxTinde)
						call Res%pointer(sRp,i)
						call sStoreSomeData(sRp,indexend,sTp,TDim,indexstart,indexend,minmaxindex(1,:),minmaxindex(2,:))
					end if
				end do
			case (3)
				call T%pointer(dTp)
				do i=1,Res%getTotalBlock()
					if(Res%getFlag(i))then
						call IndesToaddress(alldim,inde,i)
						minmaxindex=Res%NonSymIndex(inde)
						indexend=Res%getBlockdim(inde)
						!!inoutT(indexstart:indexend)=inoutT(minindeT:maxTinde)
						call Res%pointer(dRp,i)
						call dStoreSomeData(dRp,indexend,dTp,TDim,indexstart,indexend,minmaxindex(1,:),minmaxindex(2,:))
					end if
				end do
			case (4)
				call T%pointer(cTp)
				do i=1,Res%getTotalBlock()
					if(Res%getFlag(i))then
						call IndesToaddress(alldim,inde,i)
						minmaxindex=Res%NonSymIndex(inde)
						indexend=Res%getBlockdim(inde)
						!!inoutT(indexstart:indexend)=inoutT(minindeT:maxTinde)
						call Res%pointer(cRp,i)
						call cStoreSomeData(cRp,indexend,cTp,TDim,indexstart,indexend,minmaxindex(1,:),minmaxindex(2,:))
					end if
				end do
			case (5)
				call T%pointer(zTp)
				do i=1,Res%getTotalBlock()
					if(Res%getFlag(i))then
						call IndesToaddress(alldim,inde,i)
						minmaxindex=Res%NonSymIndex(inde)
						indexend=Res%getBlockdim(inde)
						!!inoutT(indexstart:indexend)=inoutT(minindeT:maxTinde)
						call Res%pointer(zRp,i)
						call zStoreSomeData(zRp,indexend,zTp,TDim,indexstart,indexend,minmaxindex(1,:),minmaxindex(2,:))
					end if
				end do
			case (6)
				call T%pointer(lTp)
				do i=1,Res%getTotalBlock()
					if(Res%getFlag(i))then
						call IndesToaddress(alldim,inde,i)
						minmaxindex=Res%NonSymIndex(inde)
						indexend=Res%getBlockdim(inde)
						!!inoutT(indexstart:indexend)=inoutT(minindeT:maxTinde)
						call Res%pointer(lRp,i)
						call lStoreSomeData(lRp,indexend,lTp,TDim,indexstart,indexend,minmaxindex(1,:),minmaxindex(2,:))
					end if
				end do
			case (7)
				call T%pointer(aTp)
				do i=1,Res%getTotalBlock()
					if(Res%getFlag(i))then
						call IndesToaddress(alldim,inde,i)
						minmaxindex=Res%NonSymIndex(inde)
						indexend=Res%getBlockdim(inde)
						!!inoutT(indexstart:indexend)=inoutT(minindeT:maxTinde)
						call Res%pointer(aRp,i)
						call aStoreSomeData(aRp,indexend,aTp,TDim,indexstart,indexend,minmaxindex(1,:),minmaxindex(2,:))
					end if
				end do
			case default
				call T%ClassPointer(clsTp)
				do i=1,Res%getTotalBlock()
					if(Res%getFlag(i))then
						call IndesToaddress(alldim,inde,i)
						minmaxindex=Res%NonSymIndex(inde)
						indexend=Res%getBlockdim(inde)
						!!inoutT(indexstart:indexend)=inoutT(minindeT:maxTinde)
						call Res%ClassPointer(clsRp,i)
						call ClassStoreSomeData(clsRp,indexend,clsTp,TDim,indexstart,indexend,minmaxindex(1,:),minmaxindex(2,:))
					end if
				end do
		end select
			
		if(T%getNameFlag())then
			do i=1,rank
				call Res%setName(i,T%getName(i))
			end do
		end if
		return
	end subroutine




	subroutine randomTensorElement(A,region)
		class(Tensor),intent(inout)::A
		class(*),intent(in),optional::region(:)
		real*8 ::temp_real,temp_imag,minnum,maxnum,delmum,numr,numi
		complex*8::numc
		complex*16::numz
		integer::i,linchar,k,inum
		character(len=100)::w
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		logical,pointer::lp(:)
		character(len=A%Data%DataCharacterLen),pointer::ap(:)
		if(present(region))then
			minnum=max(dselect(region(1)),dselect(region(2))) 
			maxnum=min(dselect(region(1)),dselect(region(2))) 
		else
			minnum=0d0
			maxnum=1d0
		end if
		if(.not.A%Data%getFlag())then
			call writemess('Can not random the element of the tensor,The Dataarray is empty',-1)
			call A%diminfo(.true.)
			call error_Stop
		end if
		delmum=maxnum-minnum
		select case(A%GetType())
			case (1)
				call A%Data%pointAllData(ip)
				do i=1,A%getTotalData()
					numr=randomnumber()
					numr=delmum*numr+minnum
					ip(i)=numr
				end do
			case (2)
				call A%Data%pointAllData(sp)
				do i=1,A%getTotalData()
					numr=randomnumber()
					numr=delmum*numr+minnum
					sp(i)=numr
				end do
			case (3)
				call A%Data%pointAllData(dp)
				do i=1,A%getTotalData()
					numr=randomnumber()
					numr=delmum*numr+minnum
					dp(i)=numr
				end do
			case (4)
				call A%Data%pointAllData(cp)
				do i=1,A%getTotalData()
					temp_real=randomnumber()
					temp_imag=randomnumber()
					numc=cmplx(temp_real,temp_imag,kind=4)
					numc=(numc*delmum)+cmplx(minnum,minnum,kind=4)
					cp(i)=numc
				end do
			case (5)
				call A%Data%pointAllData(zp)
				do i=1,A%getTotalData()
					temp_real=randomnumber()
					temp_imag=randomnumber()
					numz=dcmplx(temp_real,temp_imag)
					numz=(numz*delmum)+dcmplx(minnum,minnum)
					zp(i)=numz
				end do
			case (6)
				call A%Data%pointAllData(lp)
				do i=1,A%getTotalData()
					numr=randomnumber()
					if(numr.le. 0.5d0) then
						lp(i)=.true.
					else
						lp(i)=.false.
					end if
				end do
			case (7)
				if(present(region))then
					linchar=iselect(region(1))
				else
					linchar=3
				end if
				maxnum=126
				minnum=32
				delmum=maxnum-minnum
				call A%Data%pointAllData(ap)
				do i=1,A%getTotalData()
					w=''
					do k=1,linchar
						inum=(delmum+1)*randomnumber()+minnum
						w=w+char(inum)
					end do
					ap(i)=(trim(adjustl(w)))
				end do
		end select
	end subroutine


	subroutine setZeroValue(A)
		class(Tensor),intent(inout)::A
		if(.not.A%Data%getFlag()) return
		call ModifyAllValue(A%Data%ClassData,0d0,A%getTotalData())
		return
	end subroutine


	!*********************************************************
	!
	!   allocate
	!
	!*********************************************************

	subroutine NonZeroElement(totalBlockNum,dimen,indices,blockdim)
		integer,intent(inout)::totalBlockNum(:)
		type(Dimension),intent(in)::dimen
		integer,intent(inout)::indices(:),blockdim(:)
		integer::Rank,lenBlock,i,ii
		logical::goon,rulefit
		integer,pointer::maxdim(:),deg(:)
		call dimen%pointDim(maxdim)
		lenBlock=product(maxdim)
		if(size(totalBlockNum).ne.lenBlock)then
			call writemess('ERROR in NonZeroElement',-1)
			call writemess('size(totalBlockNum)='+size(totalBlockNum),-1)
			call writemess('product(dimData)='+lenBlock,-1)
			call dimen%diminfo(.true.)
			call error_stop
		end if
		rank=dimen%getRank()

		

		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			if(i.gt.lenBlock)then
				call writemess('ERROR in NonZeroElement,1',-1)
				call error_stop
			end if
			call CheckIndex(i,indices,MaxDim)
			rulefit=if_symmetry_rule(dimen,indices)
			if(rulefit)then
				call Dimen%pointDeg(Deg,1)
				totalBlockNum(i)=Deg(indices(1))
				do ii=2,Rank
					call Dimen%pointDeg(Deg,ii)
					totalBlockNum(i)=totalBlockNum(i)*Deg(indices(ii))
				end do
			else
				totalBlockNum(i)=0
			end if
			goon=index_counter(indices,maxdim)
		end do
		if(i.ne.lenBlock)then
			call writemess('ERROR in NonZeroElement,2',-1)
			call error_stop
		end if
		return
	end subroutine

	subroutine CheckIndex(ith,dim,MaxDim)
		integer,intent(in)::ith,dim(:),MaxDim(:)
		integer::checki
		if(.not.ProductTensor_output_check_flag)return
		checki=addressToIndes(dim,MaxDim)
		if(ith.ne.checki)then
			call writemess('ERROR in the checking index subroutine',-1)
			call writemess('MAxDim:',-1)
			call writemess(MaxDim,'I4',-1)
			call writemess('dim:',-1)
			call writemess(dim,'I4',-1)
			call writemess('ith='+ith)
			call error_stop
		end if
	end subroutine


	subroutine allocateTensoriType(T,dim,classtype,chalen)
		class(Tensor),intent(inout)::T
		integer,intent(in)::dim(:)
		integer,intent(in)::classtype
		integer,optional,intent(in)::chalen
		T%Dimension=dim
		call T%Data%allocate(product(dim),classtype,chalen)
		return
	end subroutine
	subroutine allocateTensorCha(T,dim,classtype,chalen)
		class(Tensor),intent(inout)::T
		integer,intent(in)::dim(:)
		character(len=*),intent(in)::classtype
		integer,optional,intent(in)::chalen
		T%Dimension=dim
		call T%Data%allocate(product(dim),classtype,chalen)
		return
	end subroutine

	subroutine allocateTensoriType2(T,dim,classtype,chalen)
		class(Tensor),intent(inout)::T
		type(Dimension),intent(in)::dim
		integer,intent(in)::classtype
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:),workingindices(:),workingblockdim(:)
		integer::totalBlock,rank
		T%Dimension=dim
		if(dim%GetSymmetryFlag())then
			rank=dim%getRank()
			call dim%pointDim(dimData)
			totalBlock=product(dimData)
			call WorkingMemory%Check()
			call WorkingMemory%allocate(1,totalBlock+rank+rank)
			call WorkingMemory%get_memory(BlockNum,totalBlock)
			call WorkingMemory%get_memory(workingindices,rank)
			call WorkingMemory%get_memory(workingblockdim,rank)
			call NonZeroElement(BlockNum,dim,workingindices,workingblockdim)
			call T%Data%allocate(BlockNum,classtype,chalen)
			call WorkingMemory%free()
		else
			call dim%pointDim(dimData)
			totalBlock=product(dimData)
			call T%Data%allocate(totalBlock,classtype,chalen)
		end if
		return
	end subroutine

	subroutine allocateTensorCha2(T,dim,classtype,chalen)
		class(Tensor),intent(inout)::T
		type(Dimension),intent(in)::dim
		character(len=*),intent(in)::classtype
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:),workingindices(:),workingblockdim(:)
		integer::totalBlock,rank
		T%Dimension=dim
		if(dim%GetSymmetryFlag())then
			call dim%pointDim(dimData)
			rank=dim%getRank()
			totalBlock=product(dimData)
			call WorkingMemory%Check()
			call WorkingMemory%allocate(1,totalBlock+rank+rank)
			call WorkingMemory%get_memory(BlockNum,totalBlock)
			call WorkingMemory%get_memory(workingindices,rank)
			call WorkingMemory%get_memory(workingblockdim,rank)
			call NonZeroElement(BlockNum,dim,workingindices,workingblockdim)
			call T%Data%allocate(BlockNum,classtype)
			call WorkingMemory%free()
		else
			call dim%pointDim(dimData)
			totalBlock=product(dimData)
			call T%Data%allocate(totalBlock,classtype)
		end if
		return
	end subroutine

	subroutine allocateTensoriQN(T,QN,classtype,chalen)
		class(Tensor),intent(inout)::T
		type(QuanNum),intent(in)::QN(:)
		integer,intent(in)::classtype
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:),workingindices(:),workingblockdim(:)
		integer::totalBlock,rank
		T%Dimension=QN
		rank=size(QN)
		call T%Dimension%pointDim(dimData)
		totalBlock=product(dimData)
		call WorkingMemory%Check()
		call WorkingMemory%allocate(1,totalBlock+rank+rank)
		call WorkingMemory%get_memory(BlockNum,totalBlock)
		call WorkingMemory%get_memory(workingindices,rank)
		call WorkingMemory%get_memory(workingblockdim,rank)
		call NonZeroElement(BlockNum,T%Dimension,workingindices,workingblockdim)
		call T%Data%allocate(BlockNum,classtype,chalen)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine allocateTensorChaQN(T,QN,classtype,chalen)
		class(Tensor),intent(inout)::T
		type(QuanNum),intent(in)::QN(:)
		character(len=*),intent(in)::classtype
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:),workingindices(:),workingblockdim(:)
		integer::totalBlock,rank
		T%Dimension=QN
		rank=size(QN)
		call T%Dimension%pointDim(dimData)
		totalBlock=product(dimData)
		call WorkingMemory%Check()
		call WorkingMemory%allocate(1,totalBlock+rank+rank)
		call WorkingMemory%get_memory(BlockNum,totalBlock)
		call WorkingMemory%get_memory(workingindices,rank)
		call WorkingMemory%get_memory(workingblockdim,rank)
		call NonZeroElement(BlockNum,T%Dimension,workingindices,workingblockdim)
		call T%Data%allocate(BlockNum,classtype,chalen)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine allocateTensoriT(T,TT)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(in)::TT
		T%Dimension=TT%Dimension
		call T%Data%allocate(TT%Data)
		return
	end subroutine

	subroutine allocateTensoriT2(T,TT,Type,chalen)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(in)::TT
		integer,intent(in)::Type
		integer,optional,intent(in)::chalen
		T%Dimension=TT%Dimension
		call T%Data%allocate(TT%Data,Type,chalen)
		return
	end subroutine

	subroutine allocateMomeryTensoriType2(T,dim,TotalLength,classtype,chalen)
		class(Tensor),intent(inout)::T
		type(Dimension),intent(in)::dim
		integer,intent(in)::classtype,TotalLength
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:)
		integer::totalBlock
		T%Dimension=dim
		call dim%pointDim(dimData)
		totalBlock=product(dimData)
		call T%Data%allocateDataArrayMomery(totalBlock,TotalLength,classtype,chalen)
		return
	end subroutine
	subroutine allocateMomeryTensoriType1(T,dim,TotalLength,classtype,chalen)
		class(Tensor),intent(inout)::T
		integer,intent(in)::dim(:)
		integer,intent(in)::classtype,TotalLength
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:)
		integer::totalBlock
		T%Dimension=dim
		totalBlock=product(dim)
		call T%Data%allocateDataArrayMomery(totalBlock,TotalLength,classtype,chalen)
		return
	end subroutine

	subroutine allocateMomeryTensoriQN(T,QN,TotalLength,classtype,chalen)
		class(Tensor),intent(inout)::T
		type(QuanNum),intent(in)::QN(:)
		integer,intent(in)::classtype,TotalLength
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:)
		integer::totalBlock
		T%Dimension=QN
		call T%Dimension%pointDim(dimData)
		totalBlock=product(dimData)
		call T%Data%allocateDataArrayMomery(totalBlock,TotalLength,classtype,chalen)
		return
	end subroutine

	subroutine allocateMomeryTensoriType3(T,dim,TotalLength,classtype,chalen)
		class(Tensor),intent(inout)::T
		type(Dimension),intent(in)::dim
		integer,intent(in)::TotalLength
		character(len=*),intent(in)::classtype
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:)
		integer::totalBlock
		T%Dimension=dim
		call dim%pointDim(dimData)
		totalBlock=product(dimData)
		call T%Data%allocateDataArrayMomery(totalBlock,TotalLength,classtype,chalen)
		return
	end subroutine
	subroutine allocateMomeryTensoriType4(T,dim,TotalLength,classtype,chalen)
		class(Tensor),intent(inout)::T
		integer,intent(in)::dim(:)
		integer,intent(in)::TotalLength
		character(len=*),intent(in)::classtype
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:)
		integer::totalBlock
		T%Dimension=dim
		totalBlock=product(dim)
		call T%Data%allocateDataArrayMomery(totalBlock,TotalLength,classtype,chalen)
		return
	end subroutine

	subroutine allocateMomeryTensoriQN2(T,QN,TotalLength,classtype,chalen)
		class(Tensor),intent(inout)::T
		type(QuanNum),intent(in)::QN(:)
		integer,intent(in)::TotalLength
		character(len=*),intent(in)::classtype
		integer,optional,intent(in)::chalen
		integer,pointer::BlockNum(:),dimData(:)
		integer::totalBlock
		T%Dimension=QN
		call T%Dimension%pointDim(dimData)
		totalBlock=product(dimData)
		call T%Data%allocateDataArrayMomery(totalBlock,TotalLength,classtype,chalen)
		return
	end subroutine

	subroutine allocateBlockMomery(T,vec)
		class(Tensor),intent(inout)::T
		integer,intent(in)::vec(:)
		integer,pointer::dim(:),Deg(:)
		integer::rank,i,BlockLength
		rank=T%getRank()
		!call WorkingMemory%Check()
		!call WorkingMemory%allocate(1,rank)
		!call WorkingMemory%get_memory(dim,rank)
		call T%pointDeg(Deg,1)
		BlockLength=Deg(vec(1))
		do i=2,rank
			call T%pointDeg(Deg,i)
			BlockLength=BlockLength*Deg(vec(i))
		end do
		call T%pointDim(dim)
		call T%Data%Set_block_momery(Dim,vec,BlockLength)
		!call WorkingMemory%free()
	end subroutine

	subroutine allocateBlockMomery2(T,TDim,vec)
		class(Tensor),intent(inout)::T
		integer,intent(in)::vec(:)
		integer,intent(in)::TDim(:)
		integer,pointer::Deg(:),dim(:)
		integer::i,BlockLength,rank,ith
		rank=T%getRank()
		call allocateCheck(TMPi1,rank)
		call allocateCheck(TMPi2,rank)
		ith=addressToIndes(vec,TDim)
		call GetBlockDimensionRoutine(T%Dimension,TMPi2(1:rank),ith,TMPi1(1:rank))
		BlockLength=product(TMPi2(1:rank))
		call T%Data%Set_block_momery(TDim,vec,BlockLength)
		return
	end subroutine

	subroutine allocateBlockMomery3(T,TDim,vec,BlockLength)
		class(Tensor),intent(inout)::T
		integer,intent(in)::vec(:),BlockLength
		integer,intent(in)::TDim(:)
		call T%Data%Set_block_momery(TDim,vec,BlockLength)
		return
	end subroutine

	subroutine allocateBlockMomery4(T,index,BlockLength)
		class(Tensor),intent(inout)::T
		integer,intent(in)::index,BlockLength
		call T%Data%Set_block_momery(index,BlockLength)
		return
	end subroutine
	subroutine allocateBlockMomery5(T,vec,BlockLength)
		class(Tensor),intent(inout)::T
		integer,intent(in)::vec(:),BlockLength
		integer,pointer::dim(:)
		call T%pointDim(dim)
		call T%Data%Set_block_momery(Dim,vec,BlockLength)
	end subroutine

	!*********************************************************
	!
	!   allocate from source 
	!
	!*********************************************************

	subroutine allocateTensorClassTpye1(T,dim,source,defineClassType)
		class(Tensor),intent(inout)::T
		integer,intent(in)::dim(:)
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		T%Dimension=dim
		call T%Data%allocateClassType(product(dim),source,defineClassType)
		return
	end subroutine

	subroutine allocateTensorClassTpye2(T,dim,source,defineClassType)
		class(Tensor),intent(inout)::T
		type(Dimension),intent(in)::dim
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		integer,pointer::BlockNum(:),dimData(:),workingindices(:),workingblockdim(:)
		integer::totalBlock,rank
		T%Dimension=dim
		if(dim%GetSymmetryFlag())then
			rank=dim%getRank()
			call dim%pointDim(dimData)
			totalBlock=product(dimData)
			call WorkingMemory%Check()
			call WorkingMemory%allocate(1,totalBlock+rank+rank)
			call WorkingMemory%get_memory(BlockNum,totalBlock)
			call WorkingMemory%get_memory(workingindices,rank)
			call WorkingMemory%get_memory(workingblockdim,rank)
			call NonZeroElement(BlockNum,dim,workingindices,workingblockdim)
			call T%Data%allocateClassType(BlockNum,source,defineClassType)
			call WorkingMemory%free()
		else
			call dim%pointDim(dimData)
			totalBlock=product(dimData)
			call T%Data%allocateClassType(totalBlock,source,defineClassType)
		end if
		return
	end subroutine

	subroutine allocateTensorClassTpyeQN(T,QN,source,defineClassType)
		class(Tensor),intent(inout)::T
		type(QuanNum),intent(in)::QN(:)
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		integer,pointer::BlockNum(:),dimData(:),workingindices(:),workingblockdim(:)
		integer::totalBlock,rank
		T%Dimension=QN
		rank=size(QN)
		call T%Dimension%pointDim(dimData)
		totalBlock=product(dimData)
		call WorkingMemory%Check()
		call WorkingMemory%allocate(1,totalBlock+rank+rank)
		call WorkingMemory%get_memory(BlockNum,totalBlock)
		call WorkingMemory%get_memory(workingindices,rank)
		call WorkingMemory%get_memory(workingblockdim,rank)
		call NonZeroElement(BlockNum,T%Dimension,workingindices,workingblockdim)
		call T%Data%allocateClassType(BlockNum,source,defineClassType)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine allocateTensorClassTpyeT(T,TT,source,defineClassType)
		class(Tensor),intent(inout)::T
		type(Tensor),intent(in)::TT
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		T%Dimension=TT%Dimension
		call T%Data%allocateClassType(TT%Data,source,defineClassType)
		return
	end subroutine

	subroutine allocateMomeryClassTypeiType2(T,dim,TotalLength,source,defineClassType)
		class(Tensor),intent(inout)::T
		type(Dimension),intent(in)::dim
		integer,intent(in)::TotalLength
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		integer,pointer::BlockNum(:),dimData(:)
		integer::totalBlock
		T%Dimension=dim
		call dim%pointDim(dimData)
		totalBlock=product(dimData)
		call T%Data%allocateMomeryClassType(totalBlock,TotalLength,source,defineClassType)
		return
	end subroutine
	subroutine allocateMomeryClassTypeiType1(T,dim,TotalLength,source,defineClassType)
		class(Tensor),intent(inout)::T
		integer,intent(in)::dim(:)
		integer,intent(in)::TotalLength
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		integer,pointer::BlockNum(:),dimData(:)
		integer::totalBlock
		T%Dimension=dim
		totalBlock=product(dim)
		call T%Data%allocateMomeryClassType(totalBlock,TotalLength,source,defineClassType)
		return
	end subroutine

	subroutine allocateMomeryClassTypeiQN(T,QN,TotalLength,source,defineClassType)
		class(Tensor),intent(inout)::T
		type(QuanNum),intent(in)::QN(:)
		integer,intent(in)::TotalLength
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		integer,pointer::BlockNum(:),dimData(:)
		integer::totalBlock
		T%Dimension=QN
		call T%Dimension%pointDim(dimData)
		totalBlock=product(dimData)
		call T%Data%allocateMomeryClassType(totalBlock,TotalLength,source,defineClassType)
		return
	end subroutine


	!*********************************************************
	!
	!   allocate all the block momery
	!
	!*********************************************************


	subroutine allocateAllBlock1(A,dimen,iclassType)
		class(Tensor),intent(inout)::A
		type(Dimension),intent(in)::dimen
		integer,intent(in)::iclassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		call dimen%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(dimen%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(dimen%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		A%Dimension=dimen
		call A%Data%allocate(BlockN,iclassType)
		return
	end subroutine
	subroutine allocateAllBlock2(A,dimen,iclassType)
		class(Tensor),intent(inout)::A
		type(Dimension),intent(in)::dimen
		character(len=*),intent(in)::iclassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		call dimen%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(dimen%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(dimen%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		A%Dimension=dimen
		call A%Data%allocate(BlockN,iclassType)
		return
	end subroutine
	subroutine allocateAllBlock3(A,QN,iclassType)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN(:)
		integer,intent(in)::iclassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		A%Dimension=QN
		call A%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(A%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(A%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		call A%Data%allocate(BlockN,iclassType)
		return
	end subroutine
	subroutine allocateAllBlock4(A,QN,iclassType)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN(:)
		character(len=*),intent(in)::iclassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		A%Dimension=QN
		call A%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(A%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(A%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		call A%Data%allocate(BlockN,iclassType)
		return
	end subroutine
	subroutine allocateAllBlock5(A,iclassType)
		class(Tensor),intent(inout)::A
		integer,intent(in)::iclassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		if(.not.A%getDimFlag())then
			call writemess('DO not set the dimension to the tensor yet',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('DO not set the symmetry dimension to the tensor yet',-1)
			call error_stop
		end if
		call A%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(A%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(A%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		call A%Data%allocate(BlockN,iclassType)
		return
	end subroutine
	subroutine allocateAllBlock6(A,iclassType)
		class(Tensor),intent(inout)::A
		character(len=*),intent(in)::iclassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		if(.not.A%getDimFlag())then
			call writemess('DO not set the dimension to the tensor yet',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('DO not set the symmetry dimension to the tensor yet',-1)
			call error_stop
		end if
		call A%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(A%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(A%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		call A%Data%allocate(BlockN,iclassType)
		return
	end subroutine
	subroutine allocateAllBlock7(A)
		class(Tensor),intent(inout)::A
		integer::iclassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		if(.not.A%getDimFlag())then
			call writemess('DO not set the dimension to the tensor yet',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('DO not set the symmetry dimension to the tensor yet',-1)
			call error_stop
		end if
		iclassType=A%getType()
		call A%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(A%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(A%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		call A%Data%allocate(BlockN,iclassType)
		return
	end subroutine

	!*********************************************************
	!
	!   allocate all the block momery of ClassType
	!
	!*********************************************************


	subroutine allocateAllBlockClassType1(A,dimen,mold,defineClassType)
		class(Tensor),intent(inout)::A
		type(Dimension),intent(in)::dimen
		class(*),intent(in)::mold
		integer,optional,intent(in)::defineClassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		call dimen%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(dimen%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(dimen%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		A%Dimension=dimen
		call A%Data%allocateClassType(BlockN,mold,defineClassType)
		return
	end subroutine
	subroutine allocateAllBlockClassType3(A,QN,mold,defineClassType)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN(:)
		class(*),intent(in)::mold
		integer,optional,intent(in)::defineClassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		A%Dimension=QN
		call A%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(A%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(A%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		call A%Data%allocateClassType(BlockN,mold,defineClassType)
		return
	end subroutine
	subroutine allocateAllBlockClassType5(A,mold,defineClassType)
		class(Tensor),intent(inout)::A
		class(*),intent(in)::mold
		integer,optional,intent(in)::defineClassType
		integer::i,totalLenth
		integer,pointer::Dim(:)
		integer,allocatable::BlockN(:),indices(:)
		logical::goon
		if(.not.A%getDimFlag())then
			call writemess('DO not set the dimension to the tensor yet',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('DO not set the symmetry dimension to the tensor yet',-1)
			call error_stop
		end if
		call A%pointDim(dim)
		allocate(BlockN(product(dim)))
		allocate(indices(A%getRank()))
		indices=1
		goon=.true.
		i=0
		do while(goon)
			i=i+1
			BlockN(i)=product(A%getBlockDim(indices))
			goon=index_counter(indices,dim)
		end do
		call A%Data%allocateClassType(BlockN,mold,defineClassType)
		return
	end subroutine


	!*********************************************************
	!
	!   Get information
	!
	!*********************************************************

	function getTotalDataAll(A)
		integer::getTotalDataAll
		class(Tensor),intent(in)::A
		getTotalDataAll=A%Data%getTotalDataAll()
		return
	end function
	function getTotalDataith(A,ith)
		integer::getTotalDataith
		class(Tensor),intent(in)::A
		integer,intent(in)::ith
		getTotalDataith=A%Data%getTotalData(ith)
		return
	end function
	function getTotalDatavec(A,vec)
		integer::getTotalDatavec
		class(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		integer,pointer::ei(:),iep2(:,:),iep3(:,:,:),iep4(:,:,:,:),iep5(:,:,:,:,:)
		integer,pointer::si(:),isp2(:,:),isp3(:,:,:),isp4(:,:,:,:),isp5(:,:,:,:,:)
		integer::index,i
		if(size(vec).ne.A%getRank())then
			call writemess('ERROR in A%getTotalData(vec)',-1)
			call writemess('size(vec)='+size(vec),-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		call A%pointDim(ei)
		do i=1,A%getRank()
			if(vec(i).gt.ei(i))then
				call writemess('ERROR in getTotalData(vec),vec > dim',-1)
				call writemess('Vec:')
				call writemess(vec,'I4')
				call writemess('Dimension:')
				call writemess(ei,'I4')
				call error_stop
			end if
		end do
		call A%Data%pointEndi(ei)
		call A%Data%pointStarti(si)
		select case(A%getRank())
			case(1)
				getTotalDatavec=ei(vec(1))-si(vec(1))+1
			case(2)
				iep2(1:A%dim(1),1:A%dim(2))=>ei
				isp2(1:A%dim(1),1:A%dim(2))=>si
				getTotalDatavec=iep2(vec(1),vec(2))-isp2(vec(1),vec(2))+1
			case(3)
				iep3(1:A%dim(1),1:A%dim(2),1:A%dim(3))=>ei
				isp3(1:A%dim(1),1:A%dim(2),1:A%dim(3))=>si
				getTotalDatavec=iep3(vec(1),vec(2),vec(3))-isp3(vec(1),vec(2),vec(3))+1
			case(4)
				iep4(1:A%dim(1),1:A%dim(2),1:A%dim(3),1:A%dim(4))=>ei
				isp4(1:A%dim(1),1:A%dim(2),1:A%dim(3),1:A%dim(4))=>si
				getTotalDatavec=iep4(vec(1),vec(2),vec(3),vec(4))-isp4(vec(1),vec(2),vec(3),vec(4))+1
			case(5)
				iep5(1:A%dim(1),1:A%dim(2),1:A%dim(3),1:A%dim(4),1:A%dim(5))=>ei
				isp5(1:A%dim(1),1:A%dim(2),1:A%dim(3),1:A%dim(4),1:A%dim(5))=>si
				getTotalDatavec=iep5(vec(1),vec(2),vec(3),vec(4),vec(5))- &
				             isp5(vec(1),vec(2),vec(3),vec(4),vec(5))+1
			case default
				index=addressToIndes(vec,A%dim())
				getTotalDatavec=ei(index)-si(index)+1
		end select
		return
	end function
	function getTotalBlock(A)
		integer::getTotalBlock
		class(Tensor),intent(in)::A
		getTotalBlock=A%Data%getTotalBlock()
		return
	end function
	function getType(A)
		integer::getType
		class(Tensor),intent(in)::A
		getType=A%Data%getType()
		return
	end function
	function getClassType(A)
		character(len=50)::getClassType
		class(Tensor),intent(in)::A
		getClassType=A%Data%getClassType()
		return
	end function

	function getFlagith(A,ith)
		logical::getFlagith
		class(Tensor),intent(in)::A
		integer,intent(in)::ith
		getFlagith=A%Data%getFlag(ith)
	end function

	function getTensorFlag(A)
		logical::getTensorFlag
		class(Tensor),intent(in)::A
		getTensorFlag=A%Data%getFlag().or.A%dimension%getDimFlag()
	end function

	function getDataFlag(A)
		logical::getDataFlag
		class(Tensor),intent(in)::A
		getDataFlag=A%Data%getFlag()
		return
	end function



	function getFlagVec(A,vec)
		logical::getFlagVec
		class(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		integer,pointer::ei(:),ip2(:,:),ip3(:,:,:),ip4(:,:,:,:),ip5(:,:,:,:,:),Dim(:)
		integer::index
		if(size(vec).ne.A%getRank())then
			call writemess('ERROR in A%getFlag(vec)',-1)
			call writemess('size(vec)='+size(vec),-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		call A%Data%pointEndi(ei)
		call A%pointDim(Dim)
		select case(A%getRank())
			case(1)
				getFlagVec=ei(vec(1)).gt.0
			case(2)
				ip2(1:dim(1),1:dim(2))=>ei
				getFlagVec=ip2(vec(1),vec(2)).gt.0
			case(3)
				ip3(1:dim(1),1:dim(2),1:dim(3))=>ei
				getFlagVec=ip3(vec(1),vec(2),vec(3)).gt.0
			case(4)
				ip4(1:dim(1),1:dim(2),1:dim(3),1:dim(4))=>ei
				getFlagVec=ip4(vec(1),vec(2),vec(3),vec(4)).gt.0
			case(5)
				ip5(1:dim(1),1:dim(2),1:dim(3),1:dim(4),1:dim(5))=>ei
				getFlagVec=ip5(vec(1),vec(2),vec(3),vec(4),vec(5)).gt.0
			case default
				index=addressToIndes(vec,dim)
				getFlagVec=A%Data%getFlag(index)
		end select
		return
	end function

	function getFlagDimVec(A,Dim,vec)
		logical::getFlagDimVec
		class(Tensor),intent(in)::A
		integer,intent(in)::vec(:),dim(:)
		integer,pointer::ei(:),ip2(:,:),ip3(:,:,:),ip4(:,:,:,:),ip5(:,:,:,:,:)
		integer::index
		if(size(vec).ne.size(Dim))then
			call writemess('ERROR in A%getFlag(Dim,vec)',-1)
			call writemess('size(vec)='+size(vec),-1)
			call writemess('size(Dim)='+size(Dim),-1)
			call error_stop
		end if
		call A%Data%pointEndi(ei)
		select case(size(dim))
			case(1)
				getFlagDimVec=ei(vec(1)).gt.0
			case(2)
				ip2(1:dim(1),1:dim(2))=>ei
				getFlagDimVec=ip2(vec(1),vec(2)).gt.0
			case(3)
				ip3(1:dim(1),1:dim(2),1:dim(3))=>ei
				getFlagDimVec=ip3(vec(1),vec(2),vec(3)).gt.0
			case(4)
				ip4(1:dim(1),1:dim(2),1:dim(3),1:dim(4))=>ei
				getFlagDimVec=ip4(vec(1),vec(2),vec(3),vec(4)).gt.0
			case(5)
				ip5(1:dim(1),1:dim(2),1:dim(3),1:dim(4),1:dim(5))=>ei
				getFlagDimVec=ip5(vec(1),vec(2),vec(3),vec(4),vec(5)).gt.0
			case default
				index=addressToIndes(vec,dim)
				getFlagDimVec=A%Data%getFlag(index)
		end select
		return
	end function


	function GetBlockDimension(A,ith)
		integer,allocatable::GetBlockDimension(:)
		class(Tensor),intent(in)::A
		integer,intent(in)::ith
		integer,pointer::index(:)
		integer::rank
		if(.not.A%getFlag())then
			call writemess('ERROR, The Tensor is empty',-1)
			call error_stop
		end if
		rank=A%getRank()
		allocate(GetBlockDimension(rank))
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank)
		call WorkingMemory%get_memory(index,rank)
		call IndesToaddress(A%dim(),index,ith)
		GetBlockDimension=A%getBlockDim(index)
		call WorkingMemory%free()
		return
	end function

	subroutine GetBlockDimensionRoutine(A,outDim,ith,workingindex)
		type(Dimension),intent(in)::A
		integer,intent(inout)::outDim(:)
		integer,intent(in)::ith
		integer,intent(inout)::workingindex(:)
		integer::rank
		if(.not.A%getDimFlag())then
			call writemess('ERROR, The dimension is empty',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			if(ith.ne.1)then
				call writemess('This is not a Symmetry tensor',-1)
				call error_stop
			end if
			outDim=A%dim()
			return
		end if
		rank=A%getRank()
		call IndesToaddress(A%dim(),workingindex,ith)
		outDim=A%getBlockDim(workingindex)
		return
	end subroutine

	function outAllName(T,w,typ)
		type(Tensor)::outAllName
		class(Tensor),intent(in)::T
		character(len=*),intent(in)::w
		character(len=*),intent(in),optional::typ
		integer::lenofdata
		character(len=characterlen),pointer::outchar(:)
		!character(len=Len_of_Name),pointer::outchar(:)
		!call outAllName%allocate([T%getRank()],'character',Len_of_Name)
		call outAllName%allocate([T%getRank()],'character')
		call outAllName%pointer(outchar)
		call T%outAllNameChar(outchar,lenofdata,w,typ)
		if(lenofdata.gt.0)then
			outAllName%dimension=[lenofdata]
			call outAllName%Data%resetTotalData(lenofdata)
		else
			call outAllName%empty()
		end if
		return
	end function


	subroutine hermitian_conjugate_Tensor_dimension1(A,ith)
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in hermitian_conjugate_Tensor_dimension1, the input tensor is not symmetry tensor',-1)
			call error_stop
		end if
		call hermitian_conjugate_dimension(A%Dimension,ith)
		return
	end subroutine
	subroutine hermitian_conjugate_Tensor_dimension2(A,cha)
		class(Tensor),intent(inout)::A
		character(len=*),intent(in)::cha
		integer::ith
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in hermitian_conjugate_Tensor_dimension2, the input tensor is not symmetry tensor',-1)
			call error_stop
		end if
		ith=A%FindOrder(cha)
		call hermitian_conjugate_dimension(A%Dimension,ith)
		return
	end subroutine
	subroutine hermitian_conjugate_Tensor_dimension3(A)
		class(Tensor),intent(inout)::A
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in hermitian_conjugate_Tensor_dimension3, the input tensor is not symmetry tensor',-1)
			call error_stop
		end if
		call hermitian_conjugate_dimension(A%Dimension)
		return
	end subroutine

	subroutine setDynamic(Da)
		class(Tensor),intent(inout)::Da
		call Da%Data%setDynamic()
		return
	end subroutine
	subroutine setType1(Da,classType)
		class(Tensor),intent(inout)::Da
		integer,intent(in)::classType
		call Da%Data%setType(classType)
		return
	end subroutine
	subroutine setType2(Da,classType)
		class(Tensor),intent(inout)::Da
		character(len=*),intent(in)::classType
		call Da%Data%setType(classType)
		return
	end subroutine
	logical function isnan(A)
		class(Tensor),intent(in)::A
		if(.not.A%Data%getFlag())then
			isnan=.false.
			return
		end if
		isnan=A%Data%isnanData()
		return
	end function
	logical function isnan0(A)
		class(Tensor),intent(in)::A
		if(.not.A%Data%getFlag())then
			isnan0=.false.
			return
		end if
		isnan0=A%Data%isnanData()
		return
	end function

	integer function getCharacterLen(A)
		class(Tensor),intent(in)::A
		getCharacterLen=A%Data%getCharacterLen()
		return
	end function
	!*********************************************************
	!
	! Overwrite
	!
	!*********************************************************



	subroutine deallocatableDimension(dimen)
		class(Tensor),intent(inout)::dimen
		call dimen%Dimension%deallocate()
		call dimen%Data%deallocate()
		return
	end subroutine

	subroutine emptyDimension(dimen)
		class(Tensor),intent(inout)::dimen
		call dimen%Dimension%empty()
		call dimen%Data%empty()
		return
	end subroutine

	subroutine writeExternlData(A,uni)
		class(Tensor),intent(in)::A
		integer,intent(in)::uni
		call A%dimension%write(uni)
		call A%Data%write(uni)
		return
	end subroutine

	subroutine writeExternlData2(A,w,uni)!no used
		class(Tensor),intent(in)::A
		character(len=*),intent(in)::w
		integer,intent(in)::uni
		call A%dimension%write(uni)
		call A%Data%write(uni)
		return
	end subroutine
	subroutine readExternalData(A,uni)
		class(Tensor),intent(inout)::A
		integer,intent(in)::uni
		call A%Dimension%read(uni)
		call A%Data%read(uni)
		return
	end subroutine
	subroutine readData(T,uni)
		class(Tensor),intent(inout)::T
		integer,intent(in)::uni
		CHARACTER(len=50)::notused
		integer::rank,TotalData,i,dim1,dim2,j
		integer,pointer::idata(:,:),iidata(:)
		real*4,pointer::sdata(:,:),ssdata(:)
		real*8,pointer::ddata(:,:),dddata(:)
		character(len=max_len_of_char_in_TData),pointer::adata(:,:),aadata(:)
		logical,pointer::ldata(:,:),lldata(:)
		if(T%getRank().gt.2)then
			call writemess('ERROR in reading data, only allow for rank<=2 Tensor',-1)
			call error_stop
		end if
		if(T%getRank().eq.1)then
			TotalData=T%getTotalData()
			select case(T%getType())
				case(1)
					call T%pointer(iidata)
					read(uni,*)(iidata(i),i=1,TotalData)
					nullify(iidata)
				case(2)
					call T%pointer(ssdata)
					read(uni,*)(ssdata(i),i=1,TotalData)
					nullify(ssdata)
				case(3)
					call T%pointer(dddata)
					read(uni,*)(dddata(i),i=1,TotalData)
					nullify(dddata)
				case(4)
					call writemess('ERROR in reading data, Tensor',-1)
					call writemess('Do not finished for complex data yet',-1)
					call writemess('You can real two real data and combine them into a complex one',-1)
					call writemess('for example:A and B are real*4 Tensor.',-1)
					call writemess(' call A%readData(unit1).',-1)
					call writemess(' call B%readData(unit2).',-1)
					call writemess(' C=cmplex(A,B).',-1)
					call error_stop
				case(5)
					call writemess('ERROR in reading data, Tensor',-1)
					call writemess('Do not finished for complex data yet',-1)
					call writemess('You can real two real data and combine them into a complex one',-1)
					call writemess('for example:A and B are real*8 Tensor.',-1)
					call writemess(' call A%readData(unit1).',-1)
					call writemess(' call B%readData(unit2).',-1)
					call writemess(' C=dcmplex(A,B).',-1)
					call error_stop
				case(6)
					call T%pointer(lldata)
					read(uni,*)(lldata(i),i=1,TotalData)
					nullify(lldata)
				case(7)
					call T%pointer(aadata)
					read(uni,*)(aadata(i),i=1,TotalData)
					nullify(aadata)
			end select
			return
		endif
		dim1=T%dim(1)
		dim2=T%dim(2)
		TotalData=T%getTotalData()
		select case(T%getType())
			case(1)
				call T%pointer(idata)
				do i=1,dim1
					read(uni,*)(idata(i,j),j=1,dim2)
				end do
				nullify(idata)
			case(2)
				call T%pointer(sdata)
				do i=1,dim1
					read(uni,*)(sdata(i,j),j=1,dim2)
				end do
				nullify(sdata)
			case(3)
				call T%pointer(ddata)
				do i=1,dim1
					read(uni,*)(ddata(i,j),j=1,dim2)
				end do
				nullify(ddata)
			case(4)
				call writemess('ERROR in reading data, Tensor',-1)
				call writemess('Do not finished for complex data yet',-1)
				call writemess('You can real two real data and combine them into a complex one',-1)
				call writemess('for example:A and B are real*4 Tensor.',-1)
				call writemess(' call A%readData(unit1).',-1)
				call writemess(' call B%readData(unit2).',-1)
				call writemess(' C=cmplex(A,B).',-1)
				call error_stop
			case(5)
				call writemess('ERROR in reading data, Tensor',-1)
				call writemess('Do not finished for complex data yet',-1)
				call writemess('You can real two real data and combine them into a complex one',-1)
				call writemess('for example:A and B are real*8 Tensor.',-1)
				call writemess(' call A%readData(unit1).',-1)
				call writemess(' call B%readData(unit2).',-1)
				call writemess(' C=dcmplex(A,B).',-1)
				call error_stop
			case(6)
				call T%pointer(ldata)
					do i=1,dim1
						read(uni,*)(ldata(i,j),j=1,dim2)
					end do
					nullify(ldata)
			case(7)
				call T%pointer(adata)
				do i=1,dim1
					read(uni,*)(adata(i,j),j=1,dim2)
				end do
				nullify(adata)
		end select
		return
	end subroutine

	subroutine allocateDimension(dimen,dimenData)
		class(Tensor),intent(inout)::dimen
		integer,intent(in)::dimenData(:)
		integer::classType
		classType=dimen%Data%getType()
		if(classType.eq.0)then
			call writemess('ERROR in Tensor%allocate(dimen)',-1)
			call writemess('DO not set class type yet',-1)
			call error_Stop
		end if
		call allocateTensoriType(dimen,dimenData,classtype)
		return
	end subroutine

	subroutine writemess_Tensor(mess,cpu_number)!overwrite writemess
		type(Tensor),intent(in)::mess
		integer,optional,intent(in)::cpu_number
		select case(mess%getType())
			case (1)
				call writemess(mess%ii(),cpu_number)
			case (2)
				call writemess(mess%si(),cpu_number)
			case (3)
				call writemess(mess%di(),cpu_number)
			case (4)
				call writemess(mess%ci(),cpu_number)
			case (5)
				call writemess(mess%zi(),cpu_number)
			case (6)
				call writemess(mess%li(),cpu_number)
			case (7)
				call writemess(mess%ai(),cpu_number)
		end select
		return
	end subroutine
	subroutine writemess_Tensor_form(mess,form,cpu_number)!overwrite writemess
		type(Tensor),intent(in)::mess
		character(len=*),intent(in)::form
		integer,optional,intent(in)::cpu_number
		select case(mess%getType())
			case (1)
				call writemess(mess%ii(),form,cpu_number)
			case (2)
				call writemess(mess%si(),form,cpu_number)
			case (3)
				call writemess(mess%di(),form,cpu_number)
			case (4)
				call writemess(mess%ci(),form,cpu_number)
			case (5)
				call writemess(mess%zi(),form,cpu_number)
			case (6)
				call writemess(mess%li(),form,cpu_number)
			case (7)
				call writemess(mess%ai(),form,cpu_number)
		end select
		return
	end subroutine


	!*********************************************************
	!
	! read from old verison tensor
	!
	!*********************************************************


	subroutine readOldVersionSymTensor(dimen,uni)
		class(Tensor),intent(inout) :: dimen
		integer,intent(in)::uni
		CHARACTER(len=50)::notused,classTypeChar
		integer::rank,TotalData,i,ith,TotalBlock
		logical::Flag
		type(Tensor),allocatable::Block(:)
		integer,allocatable::BlockNum(:)
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		logical,pointer::lp(:)
		character(len=CharacterLen),pointer::ap(:)
		read(uni,*)notused
		read(uni,*)notused
		
		if(notused.ne.'readable')then
			call writemess("error in reading SymTensor")
			call error_stop()
		end if
		read(uni,*)notused
		read(uni,*)notused,Flag
		if(.not.Flag)then
			read(uni,*)notused
			call dimen%empty()
			return
		end if
		read(uni,*)notused
		read(uni,*)classTypeChar
		call dimen%empty()
		read(uni,*) notused
		read(uni,*) rank
		read(uni,*) notused
		read(uni,*) TotalData
		read(uni,*) notused
		read(uni,*) TotalBlock
		read(uni,*) notused
		!call allocatecheck(dimen%block,TotalData)
		allocate(Block(TotalData))
		allocate(BlockNum(TotalData))
		BlockNum=0
		do i=1,TotalBlock
			read(uni,*)ith
			call readOldVersionTensor(block(ith),uni)
			BlockNum(ith)=block(ith)%getTotalData()
		end do
		call dimen%Data%allocate(BlockNum,classTypeChar)
		call readOldVersionSymDimen(dimen%Dimension,uni)
		read(uni,*) notused
		select case(dimen%getType())
			case(1)
				do i=1,TotalData
					if(dimen%getFlag(i))then
						call Dimen%pointer(ip,i)
						if(.not.block(i)%getFlag())then
							call writemess('ERROR',-1)
							call error_stop
						end if
						ip=block(i)
					end if
				end do
			case(2)
				do i=1,TotalData
					if(dimen%getFlag(i))then
						call Dimen%pointer(sp,i)
						if(.not.block(i)%getFlag())then
							call writemess('ERROR',-1)
							call error_stop
						end if
						sp=block(i)
					end if
				end do
			case(3)
				do i=1,TotalData
					if(dimen%getFlag(i))then
						call Dimen%pointer(dp,i)
						if(.not.block(i)%getFlag())then
							call writemess('ERROR',-1)
							call error_stop
						end if
						dp=block(i)
					end if
				end do
			case(4)
				do i=1,TotalData
					if(dimen%getFlag(i))then
						call Dimen%pointer(cp,i)
						if(.not.block(i)%getFlag())then
							call writemess('ERROR',-1)
							call error_stop
						end if
						cp=block(i)
					end if
				end do
			case(5)
				do i=1,TotalData
					if(dimen%getFlag(i))then
						call Dimen%pointer(zp,i)
						if(.not.block(i)%getFlag())then
							call writemess('ERROR',-1)
							call error_stop
						end if
						zp=block(i)
					end if
				end do
			case(6)
				do i=1,TotalData
					if(dimen%getFlag(i))then
						call Dimen%pointer(lp,i)
						if(.not.block(i)%getFlag())then
							call writemess('ERROR',-1)
							call error_stop
						end if
						lp=block(i)
					end if
				end do
			case(7)
				do i=1,TotalData
					if(dimen%getFlag(i))then
						call Dimen%pointer(ap,i)
						if(.not.block(i)%getFlag())then
							call writemess('ERROR',-1)
							call error_stop
						end if
						ap=block(i)
					end if
				end do
		end select
		return
	end subroutine

	subroutine readOldVersionTensor(T,uni)
		class(Tensor),intent(inout) :: T
		integer,intent(in)::uni
		CHARACTER(len=20)::classTypeChar
		CHARACTER(len=50)::notused
		integer::rank,TotalData,i
		integer,allocatable::idata(:)
		real*4,allocatable::sdata(:),sdata2(:)
		real*8,allocatable::ddata(:),ddata2(:)
		character(len=max_len_of_char_in_TData),allocatable::adata(:)
		logical,allocatable::ldata(:)
		type(Tensor)::TMP
		type(dimension)::dimen
		read(uni,*)notused
		read(uni,*)notused
		if(notused.ne.'readable')then
			call writemess("error in reading",-1)
			call error_stop()
		end if
		read(uni,*)notused
		read(uni,*)notused
		read(uni,*)classTypeChar
		call T%empty()
		if(classTypeChar.equ.'END')return
		read(uni,*) notused
		read(uni,*) rank
		read(uni,*) notused
		read(uni,*) TotalData
		call TMP%allocate([1],classTypeChar)
		select case(TMP%getType())
			case(1)
				allocate(idata(TotalData))
				read(uni,*) notused
				read(uni,*)(idata(i),i=1,TotalData)
				T=idata
			case(2)
				allocate(sdata(TotalData))
				read(uni,*) notused
				read(uni,*)(sdata(i),i=1,TotalData)
				T=sdata
			case(3)
				allocate(ddata(TotalData))
				read(uni,*) notused
				read(uni,*)(ddata(i),i=1,TotalData)
				T=ddata
			case(4)
				allocate(sdata(TotalData))
				allocate(sdata2(TotalData))
				read(uni,*) notused
				read(uni,*)(sdata(i),i=1,TotalData)
				read(uni,*) notused
				read(uni,*)(sdata2(i),i=1,TotalData)
				T=cmplx(sdata,sdata2,kind=4)
			case(5)
				allocate(ddata(TotalData))
				allocate(ddata2(TotalData))
				read(uni,*) notused
				read(uni,*)(ddata(i),i=1,TotalData)
				read(uni,*) notused
				read(uni,*)(ddata2(i),i=1,TotalData)
				T=dcmplx(ddata,ddata2)
			case(6)
				allocate(ldata(TotalData))
				read(uni,*) notused
				read(uni,*)(ldata(i),i=1,TotalData)
				T=ldata
			case(7)
				allocate(adata(TotalData))
				read(uni,*) notused
				read(uni,*)(adata(i),i=1,TotalData)
				T=adata
		end select
		call readOldVersionDimen(dimen,uni)
		T%dimension=dimen
		read(uni,*) notused
		return
	end subroutine

	subroutine readOldVersionSymDimen(dimen,uni)
		type(dimension),intent(inout)::dimen
		integer,intent(in)::uni
		character*50::notused
		integer::i
		type(QuanNum),allocatable::QN(:)
		type(dimension)::dimen0
		logical::QNflag
		call readOldVersionDimen(dimen0,uni)
		read(uni,*)notused,QNflag
		if(.not.QNflag)then
			dimen=dimen0
			return
		end if
		read(uni,*)notused
		allocate(QN(dimen0%getRank()))
		do i=1,dimen0%getRank()
			read(uni,*)notused
			call QN(i)%read(uni)
		end do
		read(uni,*)notused
		dimen=QN
		if(dimen0%getNameFlag())then
			do i=1,dimen0%getRank()
				call dimen%setName(i,dimen0%getName(i))
			end do
		end if
		return
	end subroutine

	subroutine readOldVersionDimen(dimen,uni)
		type(dimension),intent(inout)::dimen
		integer,intent(in)::uni
		logical::flag
		character*50::notused
		integer::lendimData,i,boundarysize,nameflag,Dimsize
		integer,allocatable::DimenData(:),boundary(:),dimnameint(:)
		character*50,allocatable::dimname(:)
		
		read(uni,*) notused
		read(uni,*)notused,flag
		read(uni,*)notused,lendimData
		if(lendimData.le.0)then
				read(uni,*) notused
				return
			end if
		allocate(DimenData(lendimData))
		read(uni,*)(DimenData(i),i=1,lendimData)
		if(.not.flag)then
			call writemess('DO not finsihed reading this case ',-1)
			call error_stop
			read(uni,*) notused
			read(uni,*)Dimsize
			read(uni,*) notused
			read(uni,*) boundarysize
			read(uni,*) notused
			allocate(boundary(boundarysize))
			read(uni,*)(boundary(i),i=1,boundarysize)
			read(uni,*) notused
			read(uni,*)	notused
			!dimen=Dimensioniniti(boundarysize,Dimsize,boundary,DimenData)
		else
			dimen=DimenData
		end if
		read(uni,*)notused,nameflag
		if(nameflag.eq.1)then
			allocate(dimname(lendimData))
			read(uni,*)notused
			read(uni,*)(dimname(i),i=1,lendimData)
			do i=1,lendimData
				call dimen%setName(i,dimname(i))
			end do
		end if
		if(nameflag.eq.2)then
				call writemess("cannot read int name",-1)
				call error_stop()
		end if
		read(uni,*) notused
		return
	end subroutine



	!*********************************************************
	!
	! printe info
	!
	!*********************************************************


	subroutine print_as_matrix(A_,PrintType)
		class(Tensor),intent(in)::A_
		character(len=*),optional,intent(in)::PrintType
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		logical,pointer::lp(:)
		character(len=A_%Data%DataCharacterLen),pointer::ap(:)
		class(*),pointer::cls(:)
		type(Tensor)::A
		integer::i
		character(len=50)::classType
		if(.not.A_%getFlag())then
			call writemess('The Tensor is empty',-1)
			return
		end if
		classType=A_%getClassType()
		call writemess('classType :'+(' '+classType))
		if(.not.A_%getDimFlag())then
			call writemess('DataArray:')
			select case(A_%getType())
				case(1)
					do i=1,A_%getTotalBlock()
						call A_%pointer(ip,i)
						if(associated(ip))then
							call iprint(ip,[1,A_%getTotalData(i)],PrintType)
						else
							call writemess('Empty Block')
						end if
					end do
				case (2)
					do i=1,A_%getTotalBlock()
						call A_%pointer(sp,i)
						if(associated(sp))then
							call sprint(sp,[1,A_%getTotalData(i)],PrintType)
						else
							call writemess('Empty Block')
						end if
					end do
				case (3)
					do i=1,A_%getTotalBlock()
						call A_%pointer(dp,i)
						if(associated(dp))then
							call dprint(dp,[1,A_%getTotalData(i)],PrintType)
						else
							call writemess('Empty Block')
						end if
					end do
				case (4)
					do i=1,A_%getTotalBlock()
						call A_%pointer(cp,i)
						if(associated(cp))then
							call cprint(cp,[1,A_%getTotalData(i)],PrintType)
						else
							call writemess('Empty Block')
						end if
					end do
				case (5)
					do i=1,A_%getTotalBlock()
						call A_%pointer(zp,i)
						if(associated(zp))then
							call zprint(zp,[1,A_%getTotalData(i)],PrintType)
						else
							call writemess('Empty Block')
						end if
					end do
				case (6)
					do i=1,A_%getTotalBlock()
						call A_%pointer(lp,i)
						if(associated(lp))then
							call lprint(lp,[1,A_%getTotalData(i)],PrintType)
						else
							call writemess('Empty Block')
						end if
					end do
				case (7)
					do i=1,A_%getTotalBlock()
						call A_%pointer(ap,i)
						if(associated(ap))then
							call aprint(ap,[1,A_%getTotalData(i)],PrintType)
						else
							call writemess('Empty Block')
						end if
					end do
				case default
					do i=1,A_%getTotalBlock()
						call A_%ClassPointer(cls,i)
						if(associated(cls))then
							call ClassPrint(cls,[1,A_%getTotalData(i)],PrintType)
						else
							call writemess('Empty Block')
						end if
					end do
			end select
			return
		end if
		if(A_%getTotalBlock().ne.1)then
			call A_%asTensor(A)
		else
			A=A_
		end if
		select case(A%getType())
			case(1)
				call A%pointer(ip)
				call iprint(ip,A%dim(),PrintType)
			case (2)
				call A%pointer(sp)
				call sprint(sp,A%dim(),PrintType)
			case (3)
				call A%pointer(dp)
				call dprint(dp,A%dim(),PrintType)
			case (4)
				call A%pointer(cp)
				call cprint(cp,A%dim(),PrintType)
			case (5)
				call A%pointer(zp)
				call zprint(zp,A%dim(),PrintType)
			case (6)
				call A%pointer(lp)
				call lprint(lp,A%dim(),PrintType)
			case (7)
				call A%pointer(ap)
				call aprint(ap,A%dim(),PrintType)
			case default
				call A%ClassPointer(cls)
				call ClassPrint(cls,A%dim(),PrintType)
		end select
		return
	end subroutine

	subroutine print_as_Block(A,PrintType)
		class(Tensor),intent(in)::A
		character(len=*),optional,intent(in)::PrintType
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		logical,pointer::lp(:)
		character(len=A%Data%DataCharacterLen),pointer::ap(:)
		integer::i
		character(len=50)::classType
		integer,allocatable::bdim(:)
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			return
		end if
		classType=A%getClassType()
		call writemess('classType :'+(' '+classType))
		if(.not.A%getDimFlag())then
			call print_as_matrix(A,PrintType)
			return
		end if
		allocate(bdim(A%getRank()))
		select case(A%getType())
			case(1)
				do i=1,A%getTotalBlock()
					if(A%Data%getFlag(i))then
						call A%pointer(ip,i)
						call writemess("block id:")
						call IndesToaddress(A%dim(),bdim,i)
						call writemess(bdim)
						call iprint(ip,A%getBlockDim(bdim),PrintType)
					end if
				end do
			case (2)
				do i=1,A%getTotalBlock()
					if(A%Data%getFlag(i))then
						call A%pointer(sp,i)
						call writemess("block id:")
						call IndesToaddress(A%dim(),bdim,i)
						call writemess(bdim)
						call writemess("block totalData="+size(sp))
						call sprint(sp,A%getBlockDim(bdim),PrintType)
					end if
				end do
			case (3)
				do i=1,A%getTotalBlock()
					if(A%Data%getFlag(i))then
						call A%pointer(dp,i)
						call writemess("block id:")
						call IndesToaddress(A%dim(),bdim,i)
						call writemess(bdim)
						call writemess("block totalData="+size(dp))
						call dprint(dp,A%getBlockDim(bdim),PrintType)
					end if
				end do
			case (4)
				do i=1,A%getTotalBlock()
					if(A%Data%getFlag(i))then
						call A%pointer(cp,i)
						call writemess("block id:")
						call IndesToaddress(A%dim(),bdim,i)
						call writemess(bdim)
						call cprint(cp,A%getBlockDim(bdim),PrintType)
					end if
				end do
			case (5)
				do i=1,A%getTotalBlock()
					if(A%Data%getFlag(i))then
						call A%pointer(zp,i)
						call writemess("block id:")
						call IndesToaddress(A%dim(),bdim,i)
						call writemess(bdim)
						call zprint(zp,A%getBlockDim(bdim),PrintType)
					end if
				end do
			case (6)
				do i=1,A%getTotalBlock()
					if(A%Data%getFlag(i))then
						call A%pointer(lp,i)
						call writemess("block id:")
						call IndesToaddress(A%dim(),bdim,i)
						call writemess(bdim)
						call lprint(lp,A%getBlockDim(bdim),PrintType)
					end if
				end do
			case (7)
				do i=1,A%getTotalBlock()
					if(A%Data%getFlag(i))then
						call A%pointer(ap,i)
						call writemess("block id:")
						call IndesToaddress(A%dim(),bdim,i)
						call writemess(bdim)
						call aprint(ap,A%getBlockDim(bdim),PrintType)
					end if
				end do
			case default
				call writemess('ERROR data type='+A%getType(),-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine print_info1(A,PrintType)
		class(Tensor),intent(in)::A
		character(len=*),optional,intent(in)::PrintType
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		logical,pointer::lp(:)
		character(len=A%Data%DataCharacterLen),pointer::ap(:)
		integer::length


		length=A%getTotalData()
		select case(A%getType())
			case(1)
				call A%Data%pointAllData(ip)
				call iprint(ip,[length],PrintType)
			case (2)
				call A%Data%pointAllData(sp)
				call sprint(sp,[length],PrintType)
			case (3)
				call A%Data%pointAllData(dp)
				call dprint(dp,[length],PrintType)
			case (4)
				call A%Data%pointAllData(cp)
				call cprint(cp,[length],PrintType)
			case (5)
				call A%Data%pointAllData(zp)
				call zprint(zp,[length],PrintType)
			case (6)
				call A%Data%pointAllData(lp)
				call lprint(lp,[length],PrintType)
			case (7)
				call A%Data%pointAllData(ap)
				call aprint(ap,[length],PrintType)
		end select
		call writemess('Dimension Data are')
		call A%diminfo(.true.)
		call writemess('TotalData='+A%getTotalData())
		call writemess('TotalBlock='+A%getTotalBlock())
		call writemess('Data type='+A%getClassType())
		return
	end subroutine
	subroutine print_info2(A,printlen0,PrintType)
		class(Tensor),intent(in)::A
		integer,intent(in)::printlen0
		character(len=*),optional,intent(in)::PrintType
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		logical,pointer::lp(:)
		character(len=A%Data%DataCharacterLen),pointer::ap(:)
		integer::length,printlen

		length=A%getTotalData()
		printlen=min(printlen0,length)
		if(printlen.gt.0)then
			select case(A%getType())
				case(1)
					call A%Data%pointAllData(ip)
					call iprint(ip(1:printlen),[printlen],PrintType)
				case (2)
					call A%Data%pointAllData(sp)
					call sprint(sp(1:printlen),[printlen],PrintType)
				case (3)
					call A%Data%pointAllData(dp)
					call dprint(dp(1:printlen),[printlen],PrintType)
				case (4)
					call A%Data%pointAllData(cp)
					call cprint(cp(1:printlen),[printlen],PrintType)
				case (5)
					call A%Data%pointAllData(zp)
					call zprint(zp(1:printlen),[printlen],PrintType)
				case (6)
					call A%Data%pointAllData(lp)
					call lprint(lp(1:printlen),[printlen],PrintType)
				case (7)
					call A%Data%pointAllData(ap)
					call aprint(ap(1:printlen),[printlen],PrintType)
			end select
		end if
		call writemess('TotalData='+A%getTotalData())
		call writemess('TotalBlock='+A%getTotalBlock())
		call writemess('Dimension Data are')
		call A%diminfo(.true.)
		call writemess('Data type='+A%getClassType())
		return
	end subroutine


	!*********************************************************
	!
	!  Other Tools
	!
	!*********************************************************

	subroutine get_index(QN,QNi,outindex)
		real*4,intent(in)::QN(:),QNi
		integer,intent(inout)::outindex
		integer::i
		outindex=0
		do i=1,size(QN)
			if(QNi.equ.QN(i))then
				outindex=i
				return
			end if
		end do
		return
	end subroutine

	!*********************************************************
	!
	!  code for MPI 
	!
	!*********************************************************

	subroutine send_Tensor(Ten1,Ten2,ID1,ID2,MPIcommon)
		type(Tensor),intent(in)::Ten1
		type(Tensor),intent(inout)::Ten2
		integer,intent(in)::ID1,ID2
		integer,optional,intent(in)::MPIcommon
		call send_DataArray(Ten1%Data,Ten2%Data,ID1,ID2,MPIcommon)
		call send_Dimension(Ten1%Dimension,Ten2%Dimension,ID1,ID2,MPIcommon)
		return
	end subroutine
	subroutine BCAST_Tensor(Ten1,ID,MPIcommon)
		type(Tensor),intent(inout)::Ten1
		integer,intent(in)::ID
		integer,optional,intent(in)::MPIcommon
		call BCAST_Dimension(Ten1%Dimension,ID,MPIcommon)
		call BCAST_DataArra(Ten1%Data,ID,MPIcommon)
		return
	end subroutine

	subroutine ALLREDUCE_Tensor(inTensor,outTensor,OP,MPIcommon)
		type(Tensor),target,intent(in)::inTensor
		type(Tensor),target,intent(inout)::outTensor
		integer::ierr
		integer,intent(in)::op
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),mpi_comm
		type(Tensor),pointer::p1,p2
		integer::length,classtype
		integer,pointer::idata(:),idata2(:)
		real*4,pointer::sdata(:),sdata2(:)
		real*8,pointer::ddata(:),ddata2(:)
		complex*8,pointer::cdata(:),cdata2(:)
		complex*16,pointer::zdata(:),zdata2(:)
		logical,pointer::ldata(:),ldata2(:)
		character(len=inTensor%Data%DataCharacterLen),pointer::adata(:),adata2(:)
		logical::ALLgoonFlag,goonFlag

		tag=1
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		p1=>inTensor
		p2=>outTensor
		if(associated(p1,p2))then
			call writemess('ERROR in ALLREDUCE_Tensor,input Tensor and output Tensor can not be the same one',-1)
			call error_stop
		end if
		p1=>null()
		p2=>null()
			! if the Tensor is empty
		goonFlag=inTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in ALLREDUCE_Tensor,the is no date in one or some Tensors',-1)
			call error_stop
		end if
			! if the Tensor is the same data type
		classtype=inTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in ALLREDUCE_Tensor,the Data type in the Tensors are not the sames',-1)
			call error_stop
		end if
			! if the length of the Tensor is the same
		length=inTensor%getTotalData()
		call MPI_BCAST(length,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getTotalData().ne.length)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in ALLREDUCE_Tensor,the length od the Tensor is not the same',-1)
			call error_stop
		end if
		call outTensor%empty()
		call outTensor%allocate(inTensor)

		select case(inTensor%getType())
			case(1)
				call inTensor%pointer(idata)
				call outTensor%pointer(idata2)
				call MPI_ALLREDUCE(idata,idata2,length,MPI_INTEGER,OP,mpi_comm,ierr)
			case(2)
				call inTensor%pointer(sdata)
				call outTensor%pointer(sdata2)
				call MPI_ALLREDUCE(sdata,sdata2,length,MPI_real,OP,mpi_comm,ierr)
			case(3)
				call inTensor%pointer(ddata)
				call outTensor%pointer(ddata2)
				call MPI_ALLREDUCE(ddata,ddata2,length,MPI_double_precision,OP,mpi_comm,ierr)
			case(4)
				call inTensor%pointer(cdata)
				call outTensor%pointer(cdata2)
				call MPI_ALLREDUCE(cdata,cdata2,length,MPI_complex,OP,mpi_comm,ierr)
			case(5)
				call inTensor%pointer(zdata)
				call outTensor%pointer(zdata2)
				call MPI_ALLREDUCE(zdata,zdata2,length,MPI_double_complex,OP,mpi_comm,ierr)
			case(6)
				call inTensor%pointer(ldata)
				call outTensor%pointer(ldata2)
				call MPI_ALLREDUCE(ldata,ldata2,length,MPI_logical,OP,mpi_comm,ierr)
			case(7)
				call inTensor%pointer(adata)
				call outTensor%pointer(adata2)
				call MPI_ALLREDUCE(adata,adata2,length,MPI_character,OP,mpi_comm,ierr)
		end  select
		return
	end subroutine


	subroutine REDUCE_Tensor(inTensor,outTensor,OP,root,MPIcommon)
		type(Tensor),target,intent(in)::inTensor
		type(Tensor),target,intent(inout)::outTensor
		integer,intent(in)::root
		integer::ierr
		integer,intent(in)::op
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),mpi_comm
		type(Tensor),pointer::p1,p2
		integer::length,classtype
		integer,pointer::idata(:),idata2(:)
		real*4,pointer::sdata(:),sdata2(:)
		real*8,pointer::ddata(:),ddata2(:)
		complex*8,pointer::cdata(:),cdata2(:)
		complex*16,pointer::zdata(:),zdata2(:)
		logical,pointer::ldata(:),ldata2(:)
		character(len=inTensor%Data%DataCharacterLen),pointer::adata(:),adata2(:)
		logical::ALLgoonFlag,goonFlag

		tag=1
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		p1=>inTensor
		p2=>outTensor
		if(associated(p1,p2))then
			call writemess('ERROR in REDUCE_Tensor,input Tensor and output Tensor can not be the same one',-1)
			call error_stop
		end if
		p1=>null()
		p2=>null()
			! if the Tensor is empty
		goonFlag=inTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in REDUCE_Tensor,the is no date in one or some Tensors',-1)
			call error_stop
		end if
			! if the Tensor is the same data type
		classtype=inTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in REDUCE_Tensor,the Data type in the Tensors are not the sames',-1)
			call error_stop
		end if
			! if the length of the Tensor is the same
		length=inTensor%getTotalData()
		call MPI_BCAST(length,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getTotalData().ne.length)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('ERROR in REDUCE_Tensor,the length od the Tensor is not the same',-1)
			call error_stop
		end if
		call outTensor%empty()
		call outTensor%allocate(inTensor)

		select case(inTensor%getType())
			case(1)
				call inTensor%pointer(idata)
				call outTensor%pointer(idata2)
				call MPI_REDUCE(idata,idata2,length,MPI_INTEGER,OP,root,mpi_comm,ierr)
			case(2)
				call inTensor%pointer(sdata)
				call outTensor%pointer(sdata2)
				call MPI_REDUCE(sdata,sdata2,length,MPI_real,OP,root,mpi_comm,ierr)
			case(3)
				call inTensor%pointer(ddata)
				call outTensor%pointer(ddata2)
				call MPI_REDUCE(ddata,ddata2,length,MPI_double_precision,OP,root,mpi_comm,ierr)
			case(4)
				call inTensor%pointer(cdata)
				call outTensor%pointer(cdata2)
				call MPI_REDUCE(cdata,cdata2,length,MPI_complex,OP,root,mpi_comm,ierr)
			case(5)
				call inTensor%pointer(zdata)
				call outTensor%pointer(zdata2)
				call MPI_REDUCE(zdata,zdata2,length,MPI_double_complex,OP,root,mpi_comm,ierr)
			case(6)
				call inTensor%pointer(ldata)
				call outTensor%pointer(ldata2)
				call MPI_REDUCE(ldata,ldata2,length,MPI_logical,OP,root,mpi_comm,ierr)
			case(7)
				call inTensor%pointer(adata)
				call outTensor%pointer(adata2)
				call MPI_REDUCE(adata,adata2,length,MPI_character,OP,root,mpi_comm,ierr)
		end  select
		return
	end subroutine

	subroutine MPI_SUM_Tensor1(inTensor,outTensor,MPIcommon)
		type(Tensor),intent(inout)::outTensor
		type(Tensor)::inTensor
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::classtype,proID,proNum,mpi_comm,i,tag,istatus(MPI_STATUS_SIZE)
		logical::goonFlag,ALLgoonFlag,Allempty
		character*20::typechar
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		classtype=inTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call mpi_comm_rank(mpi_comm,proID,ierr)
			call mpi_comm_size(mpi_comm,proNum,ierr )
			call writemess('The data type is not the same in every cpu when calling MPI_SUM_Tensor',-1)
			tag=1
			typechar=inTensor%getclassType()
			call writemess('The data type in CPU'+proID+'is classtype='+typechar,-1)
			do i=1,proNum-1
				if(proID.eq.i)call mpi_send(typechar,20,MPI_character,0,tag,mpi_comm,ierr)
				if(proID.eq.0)call mpi_recv(typechar,20,MPI_character,i,tag,mpi_comm,istatus,ierr)
				call writemess('The data type in CPU'+i+'is classtype='+typechar,-1)
			end do
			call error_stop()
		end if
		goonFlag=inTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			if(inTensor%getFlag())then
				goonFlag=.false.
			else
				goonFlag=.true.
			end if
			call MPI_ALLREDUCE(goonFlag,Allempty,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
			if(Allempty)then
				call outTensor%empty()
				return
			end if
			call writemess('There are empty Tensor in some cpu, MPI_SUM_Tensor')
			call error_stop
		end if
		if(classtype.gt.5)then
			call writemess('The data type in Tensor can not be sum,the data type is classType='+inTensor%getclassType(),-1)
			call error_stop()
		end if
		call outTensor%empty()
		call outTensor%allocate(inTensor)
		call MPI_SUM_DataArra(inTensor%Data,outTensor%Data,MPIcommon)
		return
	end subroutine
	subroutine MPI_SUM_Tensor2(inoutTensor,MPIcommon)
		type(Tensor),intent(inout)::inoutTensor
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		type(DataArray)::temp
		integer::classtype,proID,proNum,mpi_comm,i,tag,istatus(MPI_STATUS_SIZE)
		logical::goonFlag,ALLgoonFlag,Allempty
		character*20::typechar
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		classtype=inoutTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inoutTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call mpi_comm_rank(mpi_comm,proID,ierr)
			call mpi_comm_size(mpi_comm,proNum,ierr )
			call writemess('The data type is not the same in every cpu when calling MPI_SUM_Tensor',-1)
			tag=1
			typechar=inoutTensor%getclassType()
			call writemess('The data type in CPU'+proID+'is classtype='+typechar,-1)
			do i=1,proNum-1
				if(proID.eq.i)call mpi_send(typechar,20,MPI_character,0,tag,mpi_comm,ierr)
				if(proID.eq.0)call mpi_recv(typechar,20,MPI_character,i,tag,mpi_comm,istatus,ierr)
				call writemess('The data type in CPU'+i+'is classtype='+typechar,-1)
			end do
			call error_stop()
		end if
		
		goonFlag=inoutTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			if(inoutTensor%getFlag())then
				goonFlag=.false.
			else
				goonFlag=.true.
			end if
			call MPI_ALLREDUCE(goonFlag,Allempty,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
			if(Allempty)then
				call inoutTensor%empty()
				return
			end if
			call writemess('There are empty Tensor in some cpu, MPI_SUM_Tensor')
			call error_stop
		end if
		
		if(classtype.gt.5)then
			call writemess('The data type in Tensor can not be sum,the data type is classType='+inoutTensor%getclassType(),-1)
			call error_stop()
		end if
		temp=inoutTensor%Data
		call MPI_SUM_DataArra(temp,inoutTensor%Data,MPIcommon)
		return
	end subroutine
	subroutine MPI_MAX_Tensor1(inTensor,outTensor,MPIcommon)
		type(Tensor),intent(inout)::outTensor
		type(Tensor)::inTensor
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::classtype,proID,proNum,mpi_comm,i,tag,istatus(MPI_STATUS_SIZE)
		logical::goonFlag,ALLgoonFlag
		character*20::typechar
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		classtype=inTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call mpi_comm_rank(mpi_comm,proID,ierr)
			call mpi_comm_size(mpi_comm,proNum,ierr )
			call writemess('The data type is not the same in every cpu when calling MPI_MAX_Tensor',-1)
			tag=1
			typechar=inTensor%getclassType()
			call writemess('The data type in CPU'+proID+'is classtype='+typechar,-1)
			do i=1,proNum-1
				if(proID.eq.i)call mpi_send(typechar,20,MPI_character,0,tag,mpi_comm,ierr)
				if(proID.eq.0)call mpi_recv(typechar,20,MPI_character,i,tag,mpi_comm,istatus,ierr)
				call writemess('The data type in CPU'+i+'is classtype='+typechar,-1)
			end do
			call error_stop()
		end if
		goonFlag=inTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('There are empty Tensor in some cpu, MPI_MAX_Tensor')
			call error_stop
		end if
		if(classtype.ge.4)then
			call writemess('The data type in Tensor can not Find MAX,the data type is classType='+inTensor%getclassType(),-1)
			call error_stop()
		end if
		call outTensor%empty()
		call outTensor%allocate(inTensor)
		call MPI_MAX_DataArray(inTensor%Data,outTensor%Data,MPIcommon)
		return
	end subroutine
	subroutine MPI_MAX_Tensor2(inoutTensor,MPIcommon)
		type(Tensor),intent(inout)::inoutTensor
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		type(DataArray)::temp
		integer::classtype,proID,proNum,mpi_comm,i,tag,istatus(MPI_STATUS_SIZE)
		logical::goonFlag,ALLgoonFlag
		character*20::typechar
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		classtype=inoutTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inoutTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call mpi_comm_rank(mpi_comm,proID,ierr)
			call mpi_comm_size(mpi_comm,proNum,ierr )
			call writemess('The data type is not the same in every cpu when calling MPI_MAX_Tensor',-1)
			tag=1
			typechar=inoutTensor%getclassType()
			call writemess('The data type in CPU'+proID+'is classtype='+typechar,-1)
			do i=1,proNum-1
				if(proID.eq.i)call mpi_send(typechar,20,MPI_character,0,tag,mpi_comm,ierr)
				if(proID.eq.0)call mpi_recv(typechar,20,MPI_character,i,tag,mpi_comm,istatus,ierr)
				call writemess('The data type in CPU'+i+'is classtype='+typechar,-1)
			end do
			call error_stop()
		end if
		
		goonFlag=inoutTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('There are empty Tensor in some cpu, MPI_MAX_Tensor')
			call error_stop
		end if
		
		if(classtype.ge.4)then
			call writemess('The data type in Tensor can not Find MAX,the data type is classType='+inoutTensor%getclassType(),-1)
			call error_stop()
		end if
		temp=inoutTensor%Data
		call MPI_MAX_DataArray(temp,inoutTensor%Data,MPIcommon)
		return
	end subroutine
	subroutine MPI_MIN_Tensor1(inTensor,outTensor,MPIcommon)
		type(Tensor),intent(inout)::outTensor
		type(Tensor)::inTensor
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::classtype,proID,proNum,mpi_comm,i,tag,istatus(MPI_STATUS_SIZE)
		logical::goonFlag,ALLgoonFlag
		character*20::typechar
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		classtype=inTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call mpi_comm_rank(mpi_comm,proID,ierr)
			call mpi_comm_size(mpi_comm,proNum,ierr )
			call writemess('The data type is not the same in every cpu when calling MPI_MIN_Tensor',-1)
			tag=1
			typechar=inTensor%getclassType()
			call writemess('The data type in CPU'+proID+'is classtype='+typechar,-1)
			do i=1,proNum-1
				if(proID.eq.i)call mpi_send(typechar,20,MPI_character,0,tag,mpi_comm,ierr)
				if(proID.eq.0)call mpi_recv(typechar,20,MPI_character,i,tag,mpi_comm,istatus,ierr)
				call writemess('The data type in CPU'+i+'is classtype='+typechar,-1)
			end do
			call error_stop()
		end if
		goonFlag=inTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('There are empty Tensor in some cpu, MPI_MIN_Tensor')
			call error_stop
		end if
		
		if(classtype.ge.4)then
			call writemess('The data type in Tensor can not Find MIN,the data type is classType='+inTensor%getclassType(),-1)
			call error_stop()
		end if
		call outTensor%empty()
		call outTensor%allocate(inTensor)
		call MPI_MIN_DataArray(inTensor%Data,outTensor%Data,MPIcommon)
		return
	end subroutine
	subroutine MPI_MIN_Tensor2(inoutTensor,MPIcommon)
		type(Tensor),intent(inout)::inoutTensor
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		type(DataArray)::temp
		integer::classtype,proID,proNum,mpi_comm,i,tag,istatus(MPI_STATUS_SIZE)
		logical::goonFlag,ALLgoonFlag
		character*20::typechar
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		classtype=inoutTensor%getType()
		call MPI_BCAST(classtype,1,MPI_integer,0,mpi_comm,ierr)
		goonFlag=.true.
		if(inoutTensor%getType().ne.classtype)goonFlag=.false.
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call mpi_comm_rank(mpi_comm,proID,ierr)
			call mpi_comm_size(mpi_comm,proNum,ierr )
			call writemess('The data type is not the same in every cpu when calling MPI_MIN_Tensor',-1)
			tag=1
			typechar=inoutTensor%getclassType()
			call writemess('The data type in CPU'+proID+'is classtype='+typechar,-1)
			do i=1,proNum-1
				if(proID.eq.i)call mpi_send(typechar,20,MPI_character,0,tag,mpi_comm,ierr)
				if(proID.eq.0)call mpi_recv(typechar,20,MPI_character,i,tag,mpi_comm,istatus,ierr)
				call writemess('The data type in CPU'+i+'is classtype='+typechar,-1)
			end do
			call error_stop()
		end if
		goonFlag=inoutTensor%getFlag()
		call MPI_ALLREDUCE(goonFlag,ALLgoonFlag,1,MPI_logical,MPI_LAND,mpi_comm,ierr)
		if(.not.ALLgoonFlag)then
			call writemess('There are empty Tensor in some cpu, MPI_MIN_Tensor')
			call error_stop
		end if
		if(classtype.ge.4)then
			call writemess('The data type in Tensor can not Find MIN,the data type is classType='+inoutTensor%getclassType(),-1)
			call error_stop()
		end if
		temp=inoutTensor%Data
		call MPI_MIN_DataArray(temp,inoutTensor%Data,MPIcommon)
		return
	end subroutine



end module
