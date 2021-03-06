function [methodinfo,structs,enuminfo,ThunkLibName]=ps5000aWrapMFile
%PS5000AWRAPMFILE Create structures to define interfaces found in 'ps5000aWrap'.

%This function was generated by loadlibrary.m parser version 1.1.6.37 on Tue Oct 22 10:03:45 2013
%perl options:'ps5000aWrap.i -outfile=ps5000aWrapMFile.m -thunkfile=ps5000aWrap_thunk_pcwin64.c -header=ps5000aWrap.h'
ival={cell(1,0)}; % change 0 to the actual number of functions to preallocate the data.
structs=[];enuminfo=[];fcnNum=1;
fcns=struct('name',ival,'calltype',ival,'LHS',ival,'RHS',ival,'alias',ival,'thunkname', ival);
MfilePath=fileparts(mfilename('fullpath'));
ThunkLibName=fullfile(MfilePath,'ps5000aWrap_thunk_pcwin64');
% extern short __stdcall RunBlock ( short handle , long preTriggerSamples , long postTriggerSamples , unsigned long timebase , unsigned long segmentIndex ); 
fcns.thunkname{fcnNum}='int16int16longlongulongulongThunk';fcns.name{fcnNum}='RunBlock'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16', 'long', 'long', 'ulong', 'ulong'};fcnNum=fcnNum+1;
% extern short __stdcall GetStreamingLatestValues ( short handle ); 
fcns.thunkname{fcnNum}='int16int16Thunk';fcns.name{fcnNum}='GetStreamingLatestValues'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16'};fcnNum=fcnNum+1;
% extern unsigned long __stdcall AvailableData ( short handle , unsigned long * startIndex ); 
fcns.thunkname{fcnNum}='ulongint16voidPtrThunk';fcns.name{fcnNum}='AvailableData'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='ulong'; fcns.RHS{fcnNum}={'int16', 'ulongPtr'};fcnNum=fcnNum+1;
% extern short __stdcall AutoStopped ( short handle ); 
fcns.thunkname{fcnNum}='int16int16Thunk';fcns.name{fcnNum}='AutoStopped'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16'};fcnNum=fcnNum+1;
% extern short __stdcall IsReady ( short handle ); 
fcns.thunkname{fcnNum}='int16int16Thunk';fcns.name{fcnNum}='IsReady'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16'};fcnNum=fcnNum+1;
% extern short __stdcall IsTriggerReady ( short handle , unsigned long * triggeredAt ); 
fcns.thunkname{fcnNum}='int16int16voidPtrThunk';fcns.name{fcnNum}='IsTriggerReady'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16', 'ulongPtr'};fcnNum=fcnNum+1;
% extern short  __stdcall ClearTriggerReady ( void ); 
fcns.thunkname{fcnNum}='int16voidThunk';fcns.name{fcnNum}='ClearTriggerReady'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}=[];fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall SetTriggerConditions ( short handle , int * conditionsArray , short nConditions ); 
fcns.thunkname{fcnNum}='ulongint16voidPtrint16Thunk';fcns.name{fcnNum}='SetTriggerConditions'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='ulong'; fcns.RHS{fcnNum}={'int16', 'int32Ptr', 'int16'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall SetTriggerProperties ( short handle , int * propertiesArray , short nProperties , long autoTrig ); 
fcns.thunkname{fcnNum}='ulongint16voidPtrint16longThunk';fcns.name{fcnNum}='SetTriggerProperties'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='ulong'; fcns.RHS{fcnNum}={'int16', 'int32Ptr', 'int16', 'long'};fcnNum=fcnNum+1;
% extern PICO_STATUS __stdcall SetPulseWidthQualifier ( short handle , int * pwqConditionsArray , short nConditions , int direction , unsigned long lower , unsigned long upper , int type ); 
fcns.thunkname{fcnNum}='ulongint16voidPtrint16int32ulongulongint32Thunk';fcns.name{fcnNum}='SetPulseWidthQualifier'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='ulong'; fcns.RHS{fcnNum}={'int16', 'int32Ptr', 'int16', 'int32', 'ulong', 'ulong', 'int32'};fcnNum=fcnNum+1;
% extern void __stdcall setChannelCount ( short handle , short channelCount ); 
fcns.thunkname{fcnNum}='voidint16int16Thunk';fcns.name{fcnNum}='setChannelCount'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'int16', 'int16'};fcnNum=fcnNum+1;
% extern short __stdcall setEnabledChannels ( short handle , short * enabledChannels ); 
fcns.thunkname{fcnNum}='int16int16voidPtrThunk';fcns.name{fcnNum}='setEnabledChannels'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16', 'int16Ptr'};fcnNum=fcnNum+1;
% extern short __stdcall setAppAndDriverBuffers ( short handle , short channel , short * appBuffer , short * driverBuffer , unsigned long bufferLength ); 
fcns.thunkname{fcnNum}='int16int16int16voidPtrvoidPtrulongThunk';fcns.name{fcnNum}='setAppAndDriverBuffers'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16', 'int16', 'int16Ptr', 'int16Ptr', 'ulong'};fcnNum=fcnNum+1;
% extern short __stdcall setMaxMinAppAndDriverBuffers ( short handle , short channel , short * appMaxBuffer , short * appMinBuffer , short * driverMaxBuffer , short * driverMinBuffer , unsigned long bufferLength ); 
fcns.thunkname{fcnNum}='int16int16int16voidPtrvoidPtrvoidPtrvoidPtrulongThunk';fcns.name{fcnNum}='setMaxMinAppAndDriverBuffers'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}='int16'; fcns.RHS{fcnNum}={'int16', 'int16', 'int16Ptr', 'int16Ptr', 'int16Ptr', 'int16Ptr', 'ulong'};fcnNum=fcnNum+1;
methodinfo=fcns;