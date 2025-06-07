# -*- coding: utf-8 -*-
#Timestamp:5/8/2025 9:58:09 AM
#Head - 共用头部，包含所有功能。
 
HUMMER = globals().get('Hummer')
if HUMMER is None:
    raise RuntimeError('Can\'t find Hummer in globals() ')
#from head import*
from spredo import * 
init(HUMMER)

import math
import random
import time


"""
不要修改HEAD
"""


#head block end

##########################################################################################################################
#Functions
##########################################################################################################################

#deck

def imgs_folder(folder,img):
    return str(folder+"/"+img)

def Z_height(consumable, liquid_volume):
    Bottom_Z = {
        'Deepwell Plate, 96-well, 2.2ml, V-bottom': {
            (-10000, 30): lambda x: 0.5 + x * 0,
            (30, 450): lambda x: ((x - 30) / 41.86) + 2.2 ,
            (450, 2200): lambda x: (3 * (x-450) / 179.63) + 12.3,
        },
        'PCRBioRadHSP9601':{
            (-10000,15): lambda x: 0.5 +x*0,
            (15,55):lambda x: 0.0985 * x + 1.2575,
            (55,200): lambda x: 0.0491 * x + 4.0824
        },
        'DeepwellPlateDT7350504':{
            (-10000,10): lambda x: 0.5 + x*0,
            (10, 60): lambda x: 0.0323 * x + 0.72,
            (60, 250): lambda x: 0.0194 * x + 1.3579,
            (250, 1000): lambda x: 0.019 * x + 1.4324
        },
        "SC_micro_tube_2ml_sp100":{
            (-10000,80): lambda x:-16.8 + x*0,
            (80,100): lambda x: 0.03*x-17,
            (100,2000): lambda x: 0.019*x -16.227
        },
        "SC_micro_tube_05ml":{
            (-10000, 50): lambda x: 0.2 + x * 0,
            (50, 180): lambda x: ((x - 30) / 19) + 1.3 ,
            (180, 700): lambda x: ((x-180) / 55) + 11,
        }}

    if consumable in Bottom_Z:
        for volume_range, formula in Bottom_Z[consumable].items():
            if volume_range[0] <= liquid_volume < volume_range[1]:
                return round(formula(liquid_volume), 1)

    
    return 1

def Tip_touch(consumable, Z_offset=-2):
    tip_heights = {
        'PCRBioRadHSP9601': 15,
        'DeepwellPlateDT7350504': 29,
        'SC_micro_tube_2ml_sp100': 27,
        'SC_micro_tube_05ml': 25
    }

    tip_touch_x = {
        'PCRBioRadHSP9601': [((0, 5), 1.7), ((5, 15), 2)],
        'DeepwellPlateDT7350504': [((0, 11), 2.1), ((10, 29), 2.6)],
        'SC_micro_tube_2ml_sp100': [((0, 5), 2.1), ((5, 10), 2.4), ((10, 20), 2.8), ((20, 27), 3)],
        'SC_micro_tube_05ml': [((0, 5), 2.1), ((5, 10), 2.4), ((10, 20), 2.8), ((20, 27), 3)]
    }

    if consumable not in tip_heights:
        return {'IfTipTouch': False, 'TipTouchHeight': 2, 'TipTouchOffsetOfX': 5}

    height = tip_heights[consumable] + Z_offset
    if height < 0:
        return {'IfTipTouch': False, 'TipTouchHeight': 2, 'TipTouchOffsetOfX': 5}

    for (low, high), x_offset in tip_touch_x.get(consumable, []):
        if low <= height <= high:
            return {'IfTipTouch': True, 'TipTouchHeight': height, 'TipTouchOffsetOfX': x_offset}

    return {'IfTipTouch': False, 'TipTouchHeight': 2, 'TipTouchOffsetOfX': 5}

def single_asp_multi(samples,number_channels,maximum_capacity=1):
    #maximum aspiration dispenses per cicle.
    ch_full = []  # List to store full channel groups
    rest_samples = samples  # Remaining samples to process
        
    # Calculate the number of samples that can be processed in this batch
    batch_capacity = number_channels * maximum_capacity

    while rest_samples >= batch_capacity:
        ch_full.append(maximum_capacity) # Add a full batch of channels
        rest_samples -= batch_capacity
    
    #last channles. inferring capacity
    for n in range(maximum_capacity,0,-1):
        if rest_samples % (n*number_channels)==0 and rest_samples!=0:
            ch_full.append(n)
            break

    return ch_full


class Smart_loading:
    def __init__(self, dictionary, required_tips,mode):
        self.dictionary = dictionary
        self.required_tips = required_tips
        self.module=next(iter(self.dictionary))
        self.mode=mode
        self.Row_1={0:1,1:2,2:3,3:4,4:5,5:6,6:7,7:8}
        self.tips_location=self.t_load()


    def next_column(self):
        full_columns=self.dictionary[self.module]
        for column in sorted(full_columns):
            if full_columns[column] >= self.required_tips:
                return column
    
    def t_load(self):
        column=self.next_column()
        row=self.Row_1[self.dictionary[self.module][column]-self.required_tips]
        self.dictionary[self.module][column]-=self.required_tips
        return {'Module':self.module, 'Col' : column, 'Row' : row,'Tips':self.mode}
    
    def find_empty_column(self):
        columns=self.dictionary[self.module]
        for column in sorted(columns):
            if columns[column] ==0:
                del self.dictionary[self.module][column]
                return {'Module':self.module, 'Col' : column, 'Row' : 1,'Tips':self.mode}

def smart_unloading(back=None):
    if back:
        unload_tips(back)
    else:
        well = '1'+random.choice(['A', 'B'])
        unload_tips({'Module' : 'POS7', 'Well':well})  

def generate_tips_dict(pos):
    tips_dict={}
    tips_dict[pos]={}
    
    for n in range(1,13):
        tips_dict[pos][n]=8
    return tips_dict


def merge_dicts(*dicts, **kwargs):
    result = {}
    for d in dicts:
        result.update(d)
    result.update(kwargs)
    return result


def splitting_samples(samples):
    s=samples
    list_samples=[]
    while s >= 8:
        list_samples.append(8) # Add a full batch of channels
        s -= 8
    if s>0:
        list_samples.append(s)
    return list_samples

class mix_variables:
    def __init__(self, reagent_dictionary,reagent_mix):
        self.name=reagent_mix
        self.consumable=reagent_dictionary[self.name]["Consumable"]
        self.loop_list=reagent_dictionary[self.name]["parameters"]["list_asp_disp"]
        self.liquid_class_asp=reagent_dictionary[self.name]["Liquid_class"]["Aspirate"]
        self.liquid_class_disp=reagent_dictionary[self.name]["Liquid_class"]["Dispense"]
        self.liquid_class_mix=reagent_dictionary[self.name]["Liquid_class"]["Mix"]
        self.vol_reagent_extra=round(reagent_dictionary[self.name]["parameters"]["vol"]*reagent_dictionary[self.name]["parameters"]["vol_extra"],1)
        self.vol_reagent=reagent_dictionary[self.name]["parameters"]["vol"]
        self.vol_total=round(reagent_dictionary[self.name]["parameters"]["Volume_plate"],1)
        self.height=Z_height(self.consumable,self.vol_total)
        self.location_asp=reagent_dictionary[self.name]["Location_asp"]
        self.location_disp=reagent_dictionary[self.name]["Location_disp"]
        self.original_col=reagent_dictionary[self.name]["Location_disp"]['Col']
        self.original_row=reagent_dictionary[self.name]["Location_disp"]['Row']
        self.loaded=0

    def set_volume(self,new_volume):
        self.vol_total=round(self.vol_total+1*(new_volume),1)
        self.height=Z_height(self.consumable,self.vol_total)
    
    def col(self):
        if self.loaded%8==0:
            self.location_disp['Col']+=1
    
    def row(self):
        if self.loaded%8==0:
            self.location_disp['Row']=1
        else:
            self.location_disp['Row']+=1

    def set_well(self,new_loaded):
        self.loaded+=new_loaded
        self.col()
        self.row()

    def reset_original(self):
        self.location_disp['Col']=self.original_col
        self.location_disp['Row']=self.original_row
        self.loaded=0

      
def Distributing_mix(reagent_dict,reagent_object,tips_to_load,mixing=False,loops=None,empty_back=None,load_back=None):
    load_tips(tips_to_load.tips_location)
    if mixing:
        vol_mixing = reagent_object.vol_total//2 if reagent_object.vol_total//2 <150 else 150
        z_mixing=Z_height(reagent_object.consumable,vol_mixing)
        mix(merge_dicts(reagent_object.location_asp,reagent_object.liquid_class_mix,{'BottomOffsetOfZ':Z_height(reagent_object.consumable,0),'SubMixLoopCounts':int(loops),'MixOffsetOfZInLoop':z_mixing,
                                                                                     'MixOffsetOfZAfterLoop':reagent_object.height+0.5,'MixLoopVolume': vol_mixing,'Tips':tips_to_load.mode}))
    
    for x in reagent_object.loop_list:
        aspirated_volume=round(reagent_object.vol_reagent_extra*x,1)
        reagent_object.set_volume(-aspirated_volume) 
        aspirate(merge_dicts(reagent_object.location_asp,reagent_object.liquid_class_asp,{'BottomOffsetOfZ':reagent_object.height,'AspirateVolume': aspirated_volume,'Tips':tips_to_load.mode}))

        for cycle in range(x):
            dispense(merge_dicts(reagent_object.liquid_class_disp,reagent_object.location_disp,{'DispenseVolume': reagent_object.vol_reagent,'Tips':tips_to_load.mode}))
            reagent_object.set_well(1)
        if empty_back:
            empty(merge_dicts(reagent_object.location_asp,{'BottomOffsetOfZ':reagent_object.height,'DispenseRateOfP':10,'DelySeconds':time_to_delay(2),'Tips':tips_to_load.mode,"IfTipTouch":False}))

    reagent_dict[reagent_object.name]["parameters"]["Volume_plate"]=reagent_object.vol_total
    
    if load_back:
        smart_unloading(tips_to_load.tips_location)
    else:
        smart_unloading()

def time_to_delay(seconds):
    hours = int(seconds) // 3600
    minutes = (int(seconds) % 3600) // 60
    secs = seconds % 60  # Keep the fraction part
    return "%02d:%02d:%04.1f" % (hours, minutes, secs)

def reagent_volumes_string(reagent_dict):
    parts = []
    for name, info in reagent_dict.items():
        volume = info["parameters"]["Volume_plate"]
        position = info["Location_asp"]["Module"]
        well = info["Location_asp"]["Well"]
        consumable = str(info["Consumable"])  # Ensure it's string-safe
        parts.append("> Add {}: {} µl\n\tconsumable: {}\n\tposition: {} at well {}".format(
            name, volume, consumable, position, well))
    return "\n".join(parts)

class Smart_dely:
    def __init__(self,duration_seconds):
        self.duration_seconds=duration_seconds
        self.starting_time=self.record_time()

    def record_time(self):
        t=int(time.time())
        return t

    def Wait(self):
        current_time=self.record_time()
        elapsed=current_time-self.starting_time
        remaining=self.duration_seconds-elapsed
        if remaining > 0:
            dely(remaining)
        else:
            pass
    
    def get_starting_time(self):
        return time.ctime(self.starting_time)

def library_creation(samples,PCR_mode,pop_up):
    global PCR_start, PCR_25_to_4, PCR_fragmentation,PCR_Ligation, tips250ul, PCR_plate, Plate_1_3ml, Empty_box, SC_micro_tube_05ml, SC_micro_tube_2ml


    unload_tips({'Module' : 'POS8'})

    def temp_to_6(): 
        temp_set(6) 
    temp_to_6_exe = parallel_block(temp_to_6)

    pcr_open_door()
    def start_PCR():
        pcr_open_door() 
        pcr_run_methods(method = PCR_start)
        pcr_run_methods(method = PCR_25_to_4)

    start_PCR_exe=parallel_block(start_PCR)

    sample_volume=35
    folder_pics="imgs_workflow/library_creation"
    imgs_of_the_workflow=[]
    for n in range(1,14):
        imgs_of_the_workflow.append("Slide"+str(n))
        

    a=0
    imgs={}
    for n in imgs_of_the_workflow:
        a+=1
        imgs[a]=imgs_folder(folder_pics,n)


    log(samples)
    list_samples=splitting_samples(samples)

    #description: fragmentation mix very viscous. Glycerol 15%.

    reagent_dict={
        "Fragmentation Mix":{"Location_asp":{'Module':"POS5","Well":"6A"},                        
                            "Location_disp":{'Module':"POS3","Col":1,"Row":1},
                            "Consumable":SC_micro_tube_05ml,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(SC_micro_tube_05ml),{"AspirateRateOfP":10,"PreAirVolume":10,"PostAirVolume":0,"IfDely":True,"DelySeconds":time_to_delay(6),"SecondRouteRate":38}),
                                            "Dispense":merge_dicts(Tip_touch(None),{'BottomOffsetOfZ':0.5,"DispenseRateOfP":10,"IfDely":True,"DelySeconds":time_to_delay(3),"SecondRouteRate":38,"IfEmptyForDispense":False,"EmptyForDispenseOffsetOfZ":15,"EmptyForDispenseDely":time_to_delay(0.5)}),
                                            "Mix":{'PreAirVolume': (10),'DispenseVolumeAfterSubmixLoop': (10),'MixLoopAspirateRate':10,'MixLoopDispenseRate':10,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(4),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":round((samples*15)*1.10,0),"vol":14.5,"list_asp_disp":single_asp_multi(samples,1,6),"vol_extra":1.07}},
        "Magnetic Beads":{"Location_asp":{'Module':"POS5","Well":"5A"},
                            "Location_disp":{'Module':"POS6","Col":1,"Row":1},
                            "Consumable":SC_micro_tube_2ml,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(SC_micro_tube_2ml),{"AspirateRateOfP":50,"PreAirVolume":10,"PostAirVolume":0,"IfDely":True,"DelySeconds":time_to_delay(1),"SecondRouteRate":80}),
                                            "Dispense":merge_dicts(Tip_touch(None),{'BottomOffsetOfZ':1,"DispenseRateOfP":50,"IfDely":True,"DelySeconds":time_to_delay(2),"SecondRouteRate":80,"IfEmptyForDispense":False,"EmptyForDispenseOffsetOfZ":15,"EmptyForDispenseDely":time_to_delay(0.5)}),
                                            "Mix":{'PreAirVolume': (15),'DispenseVolumeAfterSubmixLoop': (15),'MixLoopAspirateRate':200,'MixLoopDispenseRate':200,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(5),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":round((samples*60)*1.10,0),"vol":60,"list_asp_disp":single_asp_multi(samples,1,2),"vol_extra":1.04}},
        "Ligation Mix":{"Location_asp":{'Module':"POS5","Well":"6B"},
                            "Location_disp":{'Module':"POS3","Col":3,"Row":1},
                            "Consumable":SC_micro_tube_2ml,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(SC_micro_tube_2ml),{"AspirateRateOfP":20,"PreAirVolume":10,"PostAirVolume":0,"IfDely":True,"DelySeconds":time_to_delay(10),"SecondRouteRate":80}),
                                            "Dispense":merge_dicts(Tip_touch(None),{'BottomOffsetOfZ':0.5,"DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(5),"SecondRouteRate":80,"IfEmptyForDispense":False,"EmptyForDispenseOffsetOfZ":15,"EmptyForDispenseDely":time_to_delay(0.5)}),
                                            "Mix":{'PreAirVolume': (15),'DispenseVolumeAfterSubmixLoop': (15),'MixLoopAspirateRate':50,'MixLoopDispenseRate':50,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(10),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":round((samples*49)*1.07,0),"vol":49,"list_asp_disp":single_asp_multi(samples,1,2),"vol_extra":1.06}},
        "Elution Buffer":{"Location_asp":{'Module':"POS5","Well":"5B"},
                            "Location_disp":{'Module':"POS6","Col":9,"Row":1},
                            "Consumable":SC_micro_tube_2ml,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(None),{"AspirateRateOfP":40,"PreAirVolume":0,"PostAirVolume":0,"IfDely":True,"DelySeconds":time_to_delay(1),"SecondRouteRate":80}),
                                            "Dispense":merge_dicts(Tip_touch(None),{'BottomOffsetOfZ':0.5,"DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(1),"SecondRouteRate":80,"IfEmptyForDispense":False,"EmptyForDispenseOffsetOfZ":Z_height(Plate_1_3ml,22),"EmptyForDispenseDely":time_to_delay(0.5)}),
                                            "Mix":{'PreAirVolume': (5),'DispenseVolumeAfterSubmixLoop': (5),'MixLoopAspirateRate':100,'MixLoopDispenseRate':100,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(2),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":round((samples*32)*1.10,0),"vol":23,"list_asp_disp":single_asp_multi(samples,1,6),"vol_extra":1.06}}
                }
    if pop_up:
        imgreplace(imgs[1])
        dialog("> Prepare fresh 80% Ethanol\n> 500ul per sample\n> Place in position 6.")
    
    fea_buffer=str(round((samples*5)*1.10,0))
    fea_enzyme=str(round((samples*10)*1.10,0))
    lig_buffer=str(round((samples*25)*1.07,0))
    lig_ligase=str(round((samples*5)*1.07,0))
    lig_DNBseq=str(round((samples*5)*1.07,0))
    lig_water=str(round((samples*15)*1.07,0))
    beads_vol=str(reagent_dict['Magnetic Beads']['parameters']['Volume_plate'])
    elution_vol=str(reagent_dict['Elution Buffer']['parameters']['Volume_plate'])
    
    
    
    reagent_text="> Fragmentation Mix:\n\tFEA Buffer MGI:"+fea_buffer+"µl\n\tFEA Enzyme MGI: "+fea_enzyme+"µl\n\tConsumable:0.5ml SC Tube\n\tPosition:POS5-6A\nLigation Mix:\n\tLigation Buffer:"+lig_buffer+"µl\n\tDNA Ligase: "+lig_ligase+"µl\n\tAdapter DNBSEQ: "+lig_DNBseq+"µl\n\tWater Nuclease Free:"+lig_water+"µl\n\tConsumable: 2ml SC Tube\n\tPosition: POS5-6B\nMagnetic Beads: "+beads_vol+"µl\n\t Consumable: 2ml SC Tube\n\tPosition: POS5-5A\nElution Buffer: "+elution_vol+"µl\n\tConsumable: 2ml SC Tube\n\tPosition: POS5-5B\n\n Press Continue when reagents ready and placed"

    if pop_up:
        imgreplace(imgs[2])
        dialog(reagent_text)

    reagent_dict.update({"ethanol":{"Location_asp":{'Module':"POS6","Col":11,"Row":1},
                            "Location_disp":{'Module':"POS6","Col":1,"Row":1},
                            "Consumable":Plate_1_3ml,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(Plate_1_3ml),{"AspirateRateOfP":50,"PreAirVolume":20,"PostAirVolume":5,"IfDely":True,"DelySeconds":time_to_delay(1),"SecondRouteRate":80}),
                                            "Dispense":merge_dicts(Tip_touch(Plate_1_3ml,0),{'BottomOffsetOfZ':33,"DispenseRateOfP":150,"IfDely":True,"DelySeconds":time_to_delay(3),"SecondRouteRate":80,"IfEmptyForDispense":True,"EmptyForDispenseOffsetOfZ":33,"EmptyForDispenseDely":time_to_delay(1)}),
                                            "Mix":{'PreAirVolume': (15),'DispenseVolumeAfterSubmixLoop': (15),'MixLoopAspirateRate':100,'MixLoopDispenseRate':100,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(6),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":500,"vol":150,"list_asp_disp":[1],"vol_extra":1}}})


    tips_dict=generate_tips_dict('POS2')


    if pop_up:
        imgreplace(imgs[3])
        dialog("> In position 1:\n\t> Place the samples\n\t> Place empty PCR strips for purified adapter-ligated libraries")


    imgreplace(imgs[4 if samples <16 else 5])
    h="1 tip box" if samples <16 else "2 tip boxes"
    if pop_up:
        dialog("> Consumables:\n\t> Place "+h+"\n\t> Place 96-well Biorad PCR plate in POS3\nPress continue. The run will start")
    else:
        dialog("> Consumables:\n\t> Place "+h+"\n\t> Place 96-well Biorad PCR plate in POS3\n> check reagent mixes, and samples are placed on the deck.\nPress continue. The run will start")


    start_PCR_exe.Wait()
    
    #Fragmentation
    ##########################################################################################################################
    report(phase = 'Fragmentation', step = 'Reaction Mixture distribution（1/3）')
    ##########################################################################################################################
    imgreplace(imgs[6])
    pcr_open_door()
    fragmentation_mix=mix_variables(reagent_dict,"Fragmentation Mix")
    tip_frag=Smart_loading(tips_dict,1,1)
    Distributing_fragmentation_mix=Distributing_mix(reagent_dict,fragmentation_mix,tip_frag,mixing=True,loops=3,empty_back=True,load_back=None)
    fragmentation_mix.reset_original()

    imgreplace(imgs[7])

    report(phase = 'Fragmentation', step = 'Transfering Samples（2/3）')
    ##########################################################################################################################
    samples_column=5

    for s in list_samples:
        tips_samples=Smart_loading(tips_dict,s,1)
        load_tips(tips_samples.tips_location)
        aspirate({"Module":"POS1-1","TipsType":0,"Tips":tips_samples.mode,"Col":samples_column,"Row":1,
                "AspirateVolume":sample_volume,"BottomOffsetOfZ":0.5,"AspirateRateOfP":20,"PreAirVolume":5,"PostAirVolume":10,"IfDely":True,"DelySeconds":time_to_delay(2),"IfTipTouch":False,"TipTouchHeight":2,"TipTouchOffsetOfX":5,"SecondRouteRate":38})
        empty(merge_dicts(fragmentation_mix.location_disp,{"Tips":int(8),"BottomOffsetOfZ":Z_height(PCR_plate,fragmentation_mix.vol_reagent),"DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.5),"IfTipTouch":False,"TipTouchHeight":2,"TipTouchOffsetOfX":5,"SecondRouteRate":38}))
        mix(merge_dicts(fragmentation_mix.location_disp,Tip_touch(PCR_plate),{"Tips":int(8),"SubMixLoopCounts":11,"MixLoopVolume":fragmentation_mix.vol_reagent,
                                                        "PreAirVolume":5,"BottomOffsetOfZ":0.5,"MixOffsetOfZInLoop":Z_height(PCR_plate,sample_volume),"MixOffsetOfZAfterLoop":Z_height(PCR_plate,sample_volume+fragmentation_mix.vol_reagent),"DispenseVolumeAfterSubmixLoop":5,
                                                        "MixLoopAspirateRate":80,"MixLoopDispenseRate":30,"SubMixLoopCompletedDely":time_to_delay(2),"DispenseRateAfterSubmixLoop":20,"MixLoopAspirateDely":time_to_delay(0.5),"MixLoopDispenseDely":time_to_delay(0.5),
                                                        "SecondRouteRate":38,"DelyAfterSubmixLoopCompleted":time_to_delay(2)}))
        smart_unloading()
        fragmentation_mix.set_well(s)
        samples_column+=1
    fragmentation_mix.reset_original()

    imgreplace(imgs[8])

    report(phase = 'Fragmentation', step = 'Reaction（3/3）')
    ##########################################################################################################################

    pcr_close_door()
    def fragmentation_reaction():
        pcr_run_methods(method=PCR_fragmentation)

    fragmentation_reaction_exe=parallel_block(fragmentation_reaction)

    fragmentation_reaction_exe.Wait()
    pcr_open_door()


    report(phase = 'Ligation', step = 'Ligation Mix distribution (1/3)')
    ##########################################################################################################################

    tips_ligation=Smart_loading(tips_dict,1,1)
    Ligation_mix=mix_variables(reagent_dict,"Ligation Mix")

    Distributing_ligation_mix=Distributing_mix(reagent_dict,Ligation_mix,tips_ligation,mixing=True,loops=10,empty_back=True,load_back=None)
    Ligation_mix.reset_original()
    imgreplace(imgs[9])



    report(phase = 'Ligation', step = 'Transfering framentation product (2/3)')
    ##########################################################################################################################
    #updating liquid class for sample + fragmentation mix
    fragmentation_mix.liquid_class_asp.update(merge_dicts(Tip_touch(None),{'AspirateRateOfP':20,'SecondRouteRate': 80}))


    for s in list_samples:
        tips_s_fr=Smart_loading(tips_dict,s,1)
        load_tips(tips_s_fr.tips_location)
        aspirate(merge_dicts(fragmentation_mix.location_disp,fragmentation_mix.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume':sample_volume+fragmentation_mix.vol_reagent,'Tips':tips_s_fr.mode}))
        empty(merge_dicts(Ligation_mix.location_disp,{"Tips":int(8),"BottomOffsetOfZ":Z_height(PCR_plate,sample_volume+fragmentation_mix.vol_reagent),"DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.1),"IfTipTouch":False,"TipTouchHeight":2,"TipTouchOffsetOfX":5,"SecondRouteRate":80}))
        mix(merge_dicts(Ligation_mix.location_disp,Tip_touch(PCR_plate),{"Tips":int(8),"SubMixLoopCounts":15,"MixLoopVolume":sample_volume+fragmentation_mix.vol_reagent if sample_volume+fragmentation_mix.vol_reagent < 150 else 100,
                                                        "PreAirVolume":5,"BottomOffsetOfZ":0.5,"MixOffsetOfZInLoop":Z_height(PCR_plate,Ligation_mix.vol_reagent),"MixOffsetOfZAfterLoop":Z_height(PCR_plate,sample_volume+fragmentation_mix.vol_reagent+Ligation_mix.vol_reagent),"DispenseVolumeAfterSubmixLoop":5,
                                                        "MixLoopAspirateRate":100,"MixLoopDispenseRate":100,"SubMixLoopCompletedDely":time_to_delay(1),"DispenseRateAfterSubmixLoop":20,"MixLoopAspirateDely":time_to_delay(0.5),"MixLoopDispenseDely":time_to_delay(0.5),
                                                        "SecondRouteRate":80,"DelyAfterSubmixLoopCompleted":time_to_delay(4)}))
        smart_unloading()
        fragmentation_mix.set_well(s)
        Ligation_mix.set_well(s)

    Ligation_mix.reset_original()
    report(phase = 'Ligation', step = 'Reaction（3/3）')
    ##########################################################################################################################


    pcr_close_door()
    def Ligation_reaction():
        pcr_run_methods(method=PCR_Ligation)
    Ligation_reaction_exe=parallel_block(Ligation_reaction)
    
    dely(300 if PCR_mode else 10)
    #Ligation
    ##########################################################################################################################
    report(phase = 'Reagent Preparation', step = 'Beads')
    ##########################################################################################################################
    magnetic_beads=mix_variables(reagent_dict,"Magnetic Beads")

    tip_beads=Smart_loading(tips_dict,1,1)
    Distributing_magnetic_beads=Distributing_mix(reagent_dict,magnetic_beads,tip_beads,mixing=True,loops=20,empty_back=True,load_back=None)
    magnetic_beads.reset_original()

    imgreplace(imgs[10])
    report(phase = 'Ligation', step = 'Reaction（3/3）')
    Ligation_reaction_exe.Wait()
    pcr_open_door()
    pcr_stop_heating()
    temp_sleep()

    #Clean_up
    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Beads Binding (1/9)')
    ##########################################################################################################################

    Ligation_mix.liquid_class_asp.update(merge_dicts(Tip_touch(None)))
    mix_volume=sample_volume+fragmentation_mix.vol_reagent+Ligation_mix.vol_reagent if sample_volume+fragmentation_mix.vol_reagent+Ligation_mix.vol_reagent < 150 else 100

    waste_tips=[]
    for s in list_samples:
        tips_beads_binding=Smart_loading(tips_dict,s,1)
        load_tips(tips_beads_binding.tips_location)
        aspirate(merge_dicts(Ligation_mix.location_disp,Ligation_mix.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume':sample_volume+fragmentation_mix.vol_reagent+Ligation_mix.vol_reagent,'Tips':tips_beads_binding.mode}))
        empty(merge_dicts(magnetic_beads.location_disp,{"Tips":int(8),"BottomOffsetOfZ":Z_height(Plate_1_3ml,magnetic_beads.vol_reagent),"DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.1),"IfTipTouch":False,"TipTouchHeight":2,"TipTouchOffsetOfX":5,"SecondRouteRate":80}))
        mix(merge_dicts(magnetic_beads.location_disp,Tip_touch(Plate_1_3ml),{"Tips":int(8),"SubMixLoopCounts":15,"MixLoopVolume":mix_volume,
                                                        "PreAirVolume":5,"BottomOffsetOfZ":0.5,"MixOffsetOfZInLoop":Z_height(Plate_1_3ml,mix_volume),"MixOffsetOfZAfterLoop":Z_height(Plate_1_3ml,sample_volume+fragmentation_mix.vol_reagent+Ligation_mix.vol_reagent+magnetic_beads.vol_reagent),"DispenseVolumeAfterSubmixLoop":5,
                                                        "MixLoopAspirateRate":100,"MixLoopDispenseRate":100,"SubMixLoopCompletedDely":time_to_delay(1),"DispenseRateAfterSubmixLoop":20,"MixLoopAspirateDely":time_to_delay(0.5),"MixLoopDispenseDely":time_to_delay(0.5),
                                                        "SecondRouteRate":80,"DelyAfterSubmixLoopCompleted":time_to_delay(2)}))
        if samples > 2:
            n=tips_beads_binding.find_empty_column()
            smart_unloading(n)
            waste_tips.append(n)
        else:
            smart_unloading()
            t=Smart_loading(tips_dict,s,1)
            waste_tips.append(t.tips_location)
        Ligation_mix.set_well(s)
        magnetic_beads.set_well(s)

    dely(300 if PCR_mode else 10)
    magnetic_up(2.3,300 if PCR_mode else 10)
    magnetic_beads.reset_original()
    imgreplace(imgs[11])



    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Waste removal (2/9)')
    ##########################################################################################################################


    waste_volume=sample_volume+fragmentation_mix.vol_reagent+Ligation_mix.vol_reagent+magnetic_beads.vol_reagent
    waste_column=5
    waste_tips_iter=iter(waste_tips)


    for s in list_samples:
        tips_waste=next(waste_tips_iter)
        tips_waste_mode=tips_waste['Tips']
        load_tips(tips_waste)
        for n in range(2):
            aspirate(merge_dicts(magnetic_beads.location_disp,magnetic_beads.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume':100 if n==0 else waste_volume-100,'Tips':tips_waste_mode}))
            empty(merge_dicts(Tip_touch(Plate_1_3ml),{"Module":"POS6","TipsType":0,"Tips":tips_waste_mode,"Col":waste_column,"Row":1,"BottomOffsetOfZ":Z_height(Plate_1_3ml,100 if n==0 else waste_volume),
                "DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.1),"SecondRouteRate":80}))
        smart_unloading()
        magnetic_beads.set_well(s)
        waste_column+=1

    magnetic_beads.reset_original()


    waste_tips_iter=iter(waste_tips)
    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Ethanol_wash 1  (3/9)')
    ##########################################################################################################################

    ethanol=mix_variables(reagent_dict,"ethanol")
    ethanol.set_volume(-ethanol.vol_reagent_extra)

    dummy_n=0
    for s in list_samples:
        tips_etho=Smart_loading(tips_dict,s,1)
        load_tips(tips_etho.tips_location)
        ethanol.location_asp['Col']=11 if dummy_n==0 else 12
        aspirate(merge_dicts(ethanol.location_asp,ethanol.liquid_class_asp,{'BottomOffsetOfZ':ethanol.height,'AspirateVolume': round(ethanol.vol_reagent_extra,1),'Tips':tips_etho.mode}))
        dispense(merge_dicts(ethanol.liquid_class_disp,ethanol.location_disp,{'DispenseVolume': ethanol.vol_reagent,'Tips':tips_etho.mode}))
        ethanol.set_well(s)
        dummy_n+=1
        if samples > 2:
            n=next(waste_tips_iter)
            smart_unloading(n)
        else:
            smart_unloading()
            t=Smart_loading(tips_dict,s,1)
            waste_tips=[]
            waste_tips.append(t.tips_location)


    reagent_dict[ethanol.name]["parameters"]["Volume_plate"]=ethanol.vol_total
    ethanol_tips_iter=iter(waste_tips)
    ethanol.reset_original()

    dely(30 if PCR_mode else 10)

    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Ethanol_wash 1 removal (4/9)')
    ##########################################################################################################################
    waste_column=5
    for s in list_samples:
        tips_waste= next(ethanol_tips_iter)
        load_tips(tips_waste)
        aspirate(merge_dicts(ethanol.location_disp,ethanol.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume':round(ethanol.vol_reagent_extra,1),'Tips':tips_waste['Tips']}))
        empty(merge_dicts(Tip_touch(Plate_1_3ml),{"Module":"POS6","TipsType":0,"Tips":tips_waste['Tips'],"Col":waste_column,"Row":1,"BottomOffsetOfZ":Z_height(Plate_1_3ml,waste_volume+ethanol.vol_reagent),
                "DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.1),"SecondRouteRate":80}))
        smart_unloading()
        ethanol.set_well(s)
        waste_column+=1


    ethanol.reset_original()
    waste_tips_iter=iter(waste_tips)

    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Ethanol_wash 2 (5/9)')
    ##########################################################################################################################
    dummy_n=0
    ethanol.set_volume(-ethanol.vol_reagent_extra)
    for s in list_samples:
        tip=Smart_loading(tips_dict,s,1)
        load_tips(tip.tips_location)
        ethanol.location_asp['Col']=11 if dummy_n==0 else 12
        aspirate(merge_dicts(ethanol.location_asp,ethanol.liquid_class_asp,{'BottomOffsetOfZ':ethanol.height,'AspirateVolume': round(ethanol.vol_reagent_extra,1),'Tips':tip.mode}))
        dispense(merge_dicts(ethanol.liquid_class_disp,ethanol.location_disp,{'DispenseVolume': ethanol.vol_reagent,'Tips':tip.mode}))
        ethanol.set_well(s)
        dummy_n+=1
        if samples > 2:
            n=next(waste_tips_iter)
            smart_unloading(n)
        else:
            smart_unloading()
            t=Smart_loading(tips_dict,s,1)
            waste_tips=[]
            waste_tips.append(t.tips_location)
            
    reagent_dict[ethanol.name]["parameters"]["Volume_plate"]=ethanol.vol_total

    dely(30 if PCR_mode else 10)

    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Ethanol_wash 2 removal (6/9)')
    ##########################################################################################################################
    ethanol.reset_original()
    waste_column=5
    ethanol_tips_iter2=iter(waste_tips)

    
    for s in list_samples:
        tips_waste= next(ethanol_tips_iter2)
        load_tips(tips_waste)
        for n in range (2):
            aspirate(merge_dicts(ethanol.location_disp,ethanol.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume':round(ethanol.vol_reagent_extra,1) if n==0 else 20,'Tips':tips_waste['Tips']}))
            empty(merge_dicts(Tip_touch(Plate_1_3ml),{"Module":"POS6","TipsType":0,"Tips":tips_waste['Tips'],"Col":waste_column,"Row":1,"BottomOffsetOfZ":Z_height(Plate_1_3ml,waste_volume+ethanol.vol_reagent*2),
                    "DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.1),"SecondRouteRate":80}))
        smart_unloading()
        ethanol.set_well(s)
        waste_column+=1

    ethanol.reset_original()




    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Beads drying (7/9)')
    ##########################################################################################################################

    bead_drying=Smart_dely(180 if PCR_mode else 10)

    report(phase = 'Reagent Preparation', step = 'Elution')
    elution=mix_variables(reagent_dict,"Elution Buffer")
    tip_elution_during_drying=Smart_loading(tips_dict,1,1)
    Distributing_elution=Distributing_mix(reagent_dict,elution,tip_elution_during_drying,mixing=True,loops=1,empty_back=True,load_back=None)
    elution.reset_original()

    
    bead_drying.Wait()
    magnetic_down(0,0)

    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Eluting beads (8/9)')
    ##########################################################################################################################
    ethanol.location_disp

    tips_elution=[]
    dummy_n=1 if samples==16 else 2

    for s in list_samples:
        if dummy_n>0:
            tips_elu=Smart_loading(tips_dict,s,1) 
        else:
            tips_dict=generate_tips_dict('POS4')
            tips_elu=Smart_loading(tips_dict,s,1)
        tips=tips_elu.tips_location
        mode = tips_elu.mode
        load_tips(tips)
        aspirate(merge_dicts(elution.location_disp,elution.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume': elution.vol_reagent+1,'Tips':mode}))
        dispense(merge_dicts(ethanol.location_disp,elution.liquid_class_disp,{'DispenseVolume': elution.vol_reagent,'Tips':mode}))
        mix(merge_dicts(ethanol.location_disp,Tip_touch(Plate_1_3ml),elution.liquid_class_mix,{"Tips":int(8),"SubMixLoopCounts":15,"MixLoopVolume":elution.vol_reagent,
                                                        "BottomOffsetOfZ":0.5,"MixOffsetOfZInLoop":Z_height(Plate_1_3ml,elution.vol_reagent),"MixOffsetOfZAfterLoop":Z_height(Plate_1_3ml,elution.vol_reagent)+1}))
        
        ethanol.set_well(s)
        elution.set_well(s)
        dummy_n-=1
        tips_elution.append(tips)
        smart_unloading(tips)

    ethanol.reset_original()

    imgreplace(imgs[12])
    dely(180 if PCR_mode else 10)
    magnetic_up(2.3,180 if PCR_mode else 10)


    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Trasnfering cleaned ligation product (9/9)')
    ##########################################################################################################################

    tips_dict

    col=1
    tips_elution_iter=iter(tips_elution)
    elution.liquid_class_disp.update({"IfEmptyForDispense":True,"EmptyForDispenseOffsetOfZ":Z_height(PCR_plate,20)})
    for s in list_samples:
        tips=next(tips_elution_iter)
        mode = tips['Tips']
        load_tips(tips)
        aspirate(merge_dicts(ethanol.location_disp,elution.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume': 20,'Tips':mode}))
        dispense(merge_dicts(elution.liquid_class_disp,{"Module":"POS1-1","TipsType":0,"Tips":mode,"Col":col,"Row":1,'DispenseVolume': 20,'Tips':mode}))
        ethanol.set_well(s)
        col+=1
        smart_unloading()

    imgreplace(imgs[13])
    unload_tips({'Module' : 'POS8'})
    pcr_stop_heating()
    magnetic_down(0,0)
    dialog("End of the Library Construction\n> Purified adapter-ligated libraries located at POS1 — store at 4°C or at -20°C.\n> Dispose of used consumables appropriately.\n> Recycle only used tips (never the box)")

def library_amplification(samples,PCR_mode,pop_up):
    global PCR_start, PCR_25_to_4, PCR_fragmentation,PCR_Ligation, tips250ul, PCR_plate, Plate_1_3ml, Empty_box, SC_micro_tube_05ml, SC_micro_tube_2ml
    unload_tips({'Module' : 'POS8'})

    def temp_to_6(): 
        temp_set(6) 
    temp_to_6_exe = parallel_block(temp_to_6)

    pcr_open_door()
    def start_PCR():
        pcr_open_door() 
        pcr_run_methods(method = PCR_start)
        pcr_run_methods(method = PCR_25_to_4)

    start_PCR_exe=parallel_block(start_PCR)

    sample_volume=20
    folder_pics="imgs_workflow/library_amplification"
    imgs_of_the_workflow=[]
    for n in range(1,14):
        imgs_of_the_workflow.append("Slide"+str(n))
        

    a=0
    imgs={}
    for n in imgs_of_the_workflow:
        a+=1
        imgs[a]=imgs_folder(folder_pics,n)


    log(samples)
    list_samples=splitting_samples(samples)


    reagent_dict={
        "Amplification Mix":{"Location_asp":{'Module':"POS5","Well":"6A"},                        
                            "Location_disp":{'Module':"POS3","Col":1,"Row":1},
                            "Consumable":SC_micro_tube_05ml,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(SC_micro_tube_05ml),{"AspirateRateOfP":10,"PreAirVolume":10,"PostAirVolume":0,"IfDely":True,"DelySeconds":time_to_delay(6),"SecondRouteRate":38}),
                                            "Dispense":merge_dicts(Tip_touch(None),{'BottomOffsetOfZ':0.5,"DispenseRateOfP":10,"IfDely":True,"DelySeconds":time_to_delay(3),"SecondRouteRate":38,"IfEmptyForDispense":False,"EmptyForDispenseOffsetOfZ":15,"EmptyForDispenseDely":time_to_delay(0.5)}),
                                            "Mix":{'PreAirVolume': (10),'DispenseVolumeAfterSubmixLoop': (10),'MixLoopAspirateRate':10,'MixLoopDispenseRate':10,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(4),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":round((samples*25)*1.07,0),"vol":24.5,"list_asp_disp":single_asp_multi(samples,1,4),"vol_extra":1.04}},
        "Magnetic Beads":{"Location_asp":{'Module':"POS5","Well":"5A"},
                            "Location_disp":{'Module':"POS6","Col":1,"Row":1},
                            "Consumable":SC_micro_tube_2ml,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(SC_micro_tube_2ml),{"AspirateRateOfP":50,"PreAirVolume":10,"PostAirVolume":0,"IfDely":True,"DelySeconds":time_to_delay(1),"SecondRouteRate":80}),
                                            "Dispense":merge_dicts(Tip_touch(None),{'BottomOffsetOfZ':1,"DispenseRateOfP":50,"IfDely":True,"DelySeconds":time_to_delay(2),"SecondRouteRate":80,"IfEmptyForDispense":False,"EmptyForDispenseOffsetOfZ":15,"EmptyForDispenseDely":time_to_delay(0.5)}),
                                            "Mix":{'PreAirVolume': (15),'DispenseVolumeAfterSubmixLoop': (15),'MixLoopAspirateRate':200,'MixLoopDispenseRate':200,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(5),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":round((samples*45)*1.10,0),"vol":45,"list_asp_disp":single_asp_multi(samples,1,3),"vol_extra":1.04}},
        "Adapters":{"Location_asp":{'Module':"POS5-1","Col":1,"Row":1},
                            "Location_disp":{'Module':"POS3","Col":1,"Row":1},
                            "Consumable":PCR_plate,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(PCR_plate),{"AspirateRateOfP":10,"PreAirVolume":5,"PostAirVolume":2,"IfDely":True,"DelySeconds":time_to_delay(0.5),"SecondRouteRate":80}),
                                            "Dispense":merge_dicts(Tip_touch(None),{'BottomOffsetOfZ':1,"DispenseRateOfP":5,"IfDely":True,"DelySeconds":time_to_delay(0.5),"SecondRouteRate":38,"IfEmptyForDispense":True,"EmptyForDispenseOffsetOfZ":1,"EmptyForDispenseDely":time_to_delay(0.5)}),
                                            "Mix":{'PreAirVolume': (15),'DispenseVolumeAfterSubmixLoop': (15),'MixLoopAspirateRate':40,'MixLoopDispenseRate':40,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(3),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":6,"vol":5,"list_asp_disp":[1],"vol_extra":1}},
        "Elution Buffer":{"Location_asp":{'Module':"POS5","Well":"5B"},
                            "Location_disp":{'Module':"POS6","Col":9,"Row":1},
                            "Consumable":SC_micro_tube_2ml,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(None),{"AspirateRateOfP":40,"PreAirVolume":0,"PostAirVolume":0,"IfDely":True,"DelySeconds":time_to_delay(1),"SecondRouteRate":80}),
                                            "Dispense":merge_dicts(Tip_touch(None),{'BottomOffsetOfZ':0.5,"DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(1),"SecondRouteRate":80,"IfEmptyForDispense":False,"EmptyForDispenseOffsetOfZ":Z_height(Plate_1_3ml,22),"EmptyForDispenseDely":time_to_delay(0.5)}),
                                            "Mix":{'PreAirVolume': (5),'DispenseVolumeAfterSubmixLoop': (5),'MixLoopAspirateRate':100,'MixLoopDispenseRate':100,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(2),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":round((samples*42)*1.15,0),"vol":40,"list_asp_disp":single_asp_multi(samples,1,3),"vol_extra":1.04}},
        "ethanol":{"Location_asp":{'Module':"POS6","Col":11,"Row":1},
                            "Location_disp":{'Module':"POS6","Col":1,"Row":1},
                            "Consumable":Plate_1_3ml,
                            "Liquid_class":{"Aspirate":merge_dicts(Tip_touch(Plate_1_3ml),{"AspirateRateOfP":50,"PreAirVolume":20,"PostAirVolume":5,"IfDely":True,"DelySeconds":time_to_delay(1),"SecondRouteRate":80}),
                                            "Dispense":merge_dicts(Tip_touch(Plate_1_3ml,0),{'BottomOffsetOfZ':33,"DispenseRateOfP":150,"IfDely":True,"DelySeconds":time_to_delay(3),"SecondRouteRate":80,"IfEmptyForDispense":True,"EmptyForDispenseOffsetOfZ":33,"EmptyForDispenseDely":time_to_delay(1)}),
                                            "Mix":{'PreAirVolume': (15),'DispenseVolumeAfterSubmixLoop': (15),'MixLoopAspirateRate':100,'MixLoopDispenseRate':100,'DispenseRateAfterSubmixLoop':10,'SubMixLoopCompletedDely':time_to_delay(6),"SecondRouteRate": 80.0,}},
                            "parameters":{"Volume_plate":500,"vol":150,"list_asp_disp":[1],"vol_extra":1}}
                }
    tips_dict=generate_tips_dict('POS2')

    if pop_up:
        imgreplace(imgs[1])
        dialog("> Prepare fresh 80% Ethanol\n> 500ul per sample\n> Place in position 6.")
    
    amp_mix=str(reagent_dict["Amplification Mix"]['parameters']['Volume_plate'])
    beads_vol=str(reagent_dict['Magnetic Beads']['parameters']['Volume_plate'])
    elution_vol=str(reagent_dict["Elution Buffer"]['parameters']['Volume_plate'])
    
    reagent_text="> Amplification Mix:\n\tAmplification mix:"+amp_mix+"µl\n\tConsumable:0.5ml SC Tube\n\tPosition:POS5-6A\nMagnetic Beads: "+beads_vol+"µl\n\t Consumable: 2ml SC Tube\n\tPosition: POS5-5A\nElution Buffer: "+elution_vol+"µl\n\tConsumable: 2ml SC Tube\n\tPosition: POS5-5B\n>UDI Adapters: 6µl\n\tPosition: POS5-1 (and 2)\n\nPress Continue when reagents ready and placed"

    if pop_up:
        imgreplace(imgs[2])
        dialog(reagent_text)
    
    if pop_up:
        imgreplace(imgs[3])
        dialog("> In position 1:\n\t> Place the adapter-ligated libraries\n\t> Place empty PCR strips for amplified cleaned libraries")
    
    imgreplace(imgs[4 if samples <16 else 5])
    h="1 tip box" if samples <16 else "2 tip boxes"
    if pop_up:
        dialog("> Consumables:\n\t> Place "+h+"\n\t> Place 96-well Biorad PCR plate in POS3\nPress continue. The run will start")
    else:
        dialog("> Check:\n\t> Place "+h+"\n\t> Place 96-well Biorad PCR plate in POS3\n> check reagent mixes, samples and UDI-adapters are placed on the deck.\nPress continue. The run will start")
        

    start_PCR_exe.Wait()
    imgreplace(imgs[6])
    
    #Amplification
    ##########################################################################################################################
    report(phase = 'Amplification', step = 'Reaction Mixture distribution（1/4）')
    ##########################################################################################################################
  
    pcr_open_door()
    amplification_mix=mix_variables(reagent_dict,"Amplification Mix")
    tip_amp=Smart_loading(tips_dict,1,1)
    Distributing_amplification_mix=Distributing_mix(reagent_dict,amplification_mix,tip_amp,mixing=True,loops=3,empty_back=True,load_back=None)
    amplification_mix.reset_original()

    imgreplace(imgs[7])

    report(phase = 'Amplification', step = 'Transfering Samples（2/4）')
    ##########################################################################################################################
    samples_column=5

    for s in list_samples:
        tips_samples=Smart_loading(tips_dict,s,1)
        load_tips(tips_samples.tips_location)
        aspirate({"Module":"POS1-1","TipsType":0,"Tips":tips_samples.mode,"Col":samples_column,"Row":1,
                "AspirateVolume":sample_volume,"BottomOffsetOfZ":0.5,"AspirateRateOfP":20,"PreAirVolume":5,"PostAirVolume":10,"IfDely":True,"DelySeconds":time_to_delay(2),"IfTipTouch":False,"TipTouchHeight":2,"TipTouchOffsetOfX":5,"SecondRouteRate":38})
        empty(merge_dicts(amplification_mix.location_disp,{"Tips":int(8),"BottomOffsetOfZ":Z_height(PCR_plate,amplification_mix.vol_reagent),"DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.5),"IfTipTouch":False,"TipTouchHeight":2,"TipTouchOffsetOfX":5,"SecondRouteRate":38}))
        mix(merge_dicts(amplification_mix.location_disp,Tip_touch(PCR_plate),{"Tips":int(8),"SubMixLoopCounts":10,"MixLoopVolume":amplification_mix.vol_reagent,
                                                        "PreAirVolume":5,"BottomOffsetOfZ":0.5,"MixOffsetOfZInLoop":Z_height(PCR_plate,sample_volume),"MixOffsetOfZAfterLoop":Z_height(PCR_plate,sample_volume+amplification_mix.vol_reagent),"DispenseVolumeAfterSubmixLoop":5,
                                                        "MixLoopAspirateRate":80,"MixLoopDispenseRate":30,"SubMixLoopCompletedDely":time_to_delay(2),"DispenseRateAfterSubmixLoop":20,"MixLoopAspirateDely":time_to_delay(0.5),"MixLoopDispenseDely":time_to_delay(0.5),
                                                        "SecondRouteRate":38,"DelyAfterSubmixLoopCompleted":time_to_delay(2)}))
        smart_unloading()
        amplification_mix.set_well(s)
        samples_column+=1
    amplification_mix.reset_original()
    
    
    report(phase = 'Amplification', step = 'Transfering Adapters（3/4）')
    ##########################################################################################################################
    samples_column=1
    adapters=mix_variables(reagent_dict,"Adapters")
    
    for s in list_samples:
        tips_samples=Smart_loading(tips_dict,s,1)
        load_tips(tips_samples.tips_location)
        aspirate({"Module":"POS5-1","TipsType":0,"Tips":tips_samples.mode,"Col":samples_column,"Row":1,
                "AspirateVolume":adapters.vol_reagent,"BottomOffsetOfZ":0.5,"AspirateRateOfP":10,"PreAirVolume":5,"PostAirVolume":2,"IfDely":True,"DelySeconds":time_to_delay(2),"IfTipTouch":False,"TipTouchHeight":2,"TipTouchOffsetOfX":5,"SecondRouteRate":38})
        empty(merge_dicts(amplification_mix.location_disp,{"Tips":tips_samples.mode,"BottomOffsetOfZ":Z_height(PCR_plate,amplification_mix.vol_reagent),"DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.5),"IfTipTouch":False,"TipTouchHeight":2,"TipTouchOffsetOfX":5,"SecondRouteRate":38}))
        mix(merge_dicts(amplification_mix.location_disp,Tip_touch(PCR_plate),{"Tips":tips_samples.mode,"SubMixLoopCounts":10,"MixLoopVolume":amplification_mix.vol_reagent,
                                                        "PreAirVolume":5,"BottomOffsetOfZ":0.5,"MixOffsetOfZInLoop":Z_height(PCR_plate,sample_volume),
                                                        "MixOffsetOfZAfterLoop":Z_height(PCR_plate,sample_volume+amplification_mix.vol_reagent+adapters.vol_reagent),"DispenseVolumeAfterSubmixLoop":5,
                                                        "MixLoopAspirateRate":80,"MixLoopDispenseRate":30,"SubMixLoopCompletedDely":time_to_delay(2),"DispenseRateAfterSubmixLoop":20,
                                                        "MixLoopAspirateDely":time_to_delay(0.5),"MixLoopDispenseDely":time_to_delay(0.5),
                                                        "SecondRouteRate":38,"DelyAfterSubmixLoopCompleted":time_to_delay(2)}))
        smart_unloading()
        amplification_mix.set_well(s)
        samples_column+=1
    amplification_mix.reset_original()

  
    imgreplace(imgs[8])

    report(phase = 'Amplification', step = 'Reaction（4/4）')
    ##########################################################################################################################

    pcr_close_door()
    def amplification_reaction():
        pcr_run_methods(method=PCR_Mix)

    amplification_reaction_exe=parallel_block(amplification_reaction)
    
    dely(700 if PCR_mode else 10)
    ##########################################################################################################################
    report(phase = 'Reagent Preparation', step = 'Beads')
    ##########################################################################################################################
    magnetic_beads=mix_variables(reagent_dict,"Magnetic Beads")

    tip_beads=Smart_loading(tips_dict,1,1)
    Distributing_magnetic_beads=Distributing_mix(reagent_dict,magnetic_beads,tip_beads,mixing=True,loops=20,empty_back=True,load_back=None)
    magnetic_beads.reset_original()
    imgreplace(imgs[9])

    amplification_reaction_exe.Wait()
    
    pcr_open_door()
    pcr_stop_heating()
    temp_sleep()
    
    
    #Clean_up
    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Beads Binding (1/9)')
    ##########################################################################################################################

    amplification_mix.liquid_class_asp.update(merge_dicts(Tip_touch(None)))
    mix_volume=sample_volume+amplification_mix.vol_reagent+adapters.vol_reagent if sample_volume+amplification_mix.vol_reagent+adapters.vol_reagent < 150 else 100



    waste_tips=[]
    for s in list_samples:
        tips_beads_binding=Smart_loading(tips_dict,s,1)
        load_tips(tips_beads_binding.tips_location)
        aspirate(merge_dicts(amplification_mix.location_disp,amplification_mix.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume':sample_volume+amplification_mix.vol_reagent+adapters.vol_reagent,'Tips':tips_beads_binding.mode}))
        empty(merge_dicts(magnetic_beads.location_disp,{"Tips":tips_beads_binding.mode,"BottomOffsetOfZ":Z_height(Plate_1_3ml,magnetic_beads.vol_reagent),"DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.1),"IfTipTouch":False,"TipTouchHeight":2,"TipTouchOffsetOfX":5,"SecondRouteRate":80}))
        mix(merge_dicts(magnetic_beads.location_disp,Tip_touch(Plate_1_3ml),{"Tips":tips_beads_binding.mode,"SubMixLoopCounts":15,"MixLoopVolume":mix_volume,
                                                        "PreAirVolume":5,"BottomOffsetOfZ":0.5,"MixOffsetOfZInLoop":Z_height(Plate_1_3ml,mix_volume),"MixOffsetOfZAfterLoop":Z_height(Plate_1_3ml,sample_volume+amplification_mix.vol_reagent+adapters.vol_reagent+magnetic_beads.vol_reagent),"DispenseVolumeAfterSubmixLoop":5,
                                                        "MixLoopAspirateRate":100,"MixLoopDispenseRate":100,"SubMixLoopCompletedDely":time_to_delay(1),"DispenseRateAfterSubmixLoop":20,"MixLoopAspirateDely":time_to_delay(0.5),"MixLoopDispenseDely":time_to_delay(0.5),
                                                        "SecondRouteRate":80,"DelyAfterSubmixLoopCompleted":time_to_delay(2)}))
        if samples > 3:
            n=tips_beads_binding.find_empty_column()
            smart_unloading(n)
            waste_tips.append(n)
        else:
            smart_unloading()
            t=Smart_loading(tips_dict,s,1)
            waste_tips.append(t.tips_location)
        amplification_mix.set_well(s)
        magnetic_beads.set_well(s)

    imgreplace(imgs[10])


    dely(300 if PCR_mode else 10)
    magnetic_up(2.3,300 if PCR_mode else 10)
    magnetic_beads.reset_original()



    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Waste removal (2/9)')
    ##########################################################################################################################


    waste_volume=mix_volume+magnetic_beads.vol_reagent
    waste_column=5
    waste_tips_iter=iter(waste_tips)


    for s in list_samples:
        tips_waste=next(waste_tips_iter)
        tips_waste_mode=tips_waste['Tips']
        load_tips(tips_waste)
        for n in range(2):
            aspirate(merge_dicts(magnetic_beads.location_disp,magnetic_beads.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume':90 if n==0 else waste_volume-70,'Tips':tips_waste_mode}))
            empty(merge_dicts(Tip_touch(Plate_1_3ml),{"Module":"POS6","TipsType":0,"Tips":tips_waste_mode,"Col":waste_column,"Row":1,"BottomOffsetOfZ":Z_height(Plate_1_3ml,100 if n==0 else waste_volume),
                "DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.1),"SecondRouteRate":80}))
        smart_unloading()
        magnetic_beads.set_well(s)
        waste_column+=1

    magnetic_beads.reset_original()


    waste_tips_iter=iter(waste_tips)
    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Ethanol_wash 1  (3/9)')
    ##########################################################################################################################

    ethanol=mix_variables(reagent_dict,"ethanol")
    ethanol.set_volume(-ethanol.vol_reagent_extra)
    
    dummy_n=0
    for s in list_samples:
        tips_etho=Smart_loading(tips_dict,s,1)
        load_tips(tips_etho.tips_location)
        ethanol.location_asp['Col']=11 if dummy_n==0 else 12
        aspirate(merge_dicts(ethanol.location_asp,ethanol.liquid_class_asp,{'BottomOffsetOfZ':ethanol.height,'AspirateVolume': round(ethanol.vol_reagent_extra,1),'Tips':tips_etho.mode}))
        dispense(merge_dicts(ethanol.liquid_class_disp,ethanol.location_disp,{'DispenseVolume': ethanol.vol_reagent,'Tips':tips_etho.mode}))
        ethanol.set_well(s)
        dummy_n+=1
        if samples > 3:
            n=next(waste_tips_iter)
            smart_unloading(n)
        else:
            smart_unloading()
            t=Smart_loading(tips_dict,s,1)
            waste_tips=[]
            waste_tips.append(t.tips_location)


    reagent_dict[ethanol.name]["parameters"]["Volume_plate"]=ethanol.vol_total
    ethanol_tips_iter=iter(waste_tips)
    ethanol.reset_original()

    dely(30 if PCR_mode else 10)

    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Ethanol_wash 1 removal (4/9)')
    ##########################################################################################################################
    waste_column=5
    for s in list_samples:
        tips_waste= next(ethanol_tips_iter)
        load_tips(tips_waste)
        aspirate(merge_dicts(ethanol.location_disp,ethanol.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume':round(ethanol.vol_reagent_extra,1),'Tips':tips_waste['Tips']}))
        empty(merge_dicts(Tip_touch(Plate_1_3ml),{"Module":"POS6","TipsType":0,"Tips":tips_waste['Tips'],"Col":waste_column,"Row":1,"BottomOffsetOfZ":Z_height(Plate_1_3ml,waste_volume+ethanol.vol_reagent),
                "DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.1),"SecondRouteRate":80}))
        smart_unloading()
        ethanol.set_well(s)
        waste_column+=1


    ethanol.reset_original()
    waste_tips_iter=iter(waste_tips)

    report(phase = 'Clean-up', step = 'Ethanol_wash 2 (5/9)')
    dummy_n=0
    ethanol.set_volume(-ethanol.vol_reagent_extra)
    for s in list_samples:
        tip=Smart_loading(tips_dict,s,1)
        load_tips(tip.tips_location)
        ethanol.location_asp['Col']=11 if dummy_n==0 else 12
        aspirate(merge_dicts(ethanol.location_asp,ethanol.liquid_class_asp,{'BottomOffsetOfZ':ethanol.height,'AspirateVolume': round(ethanol.vol_reagent_extra,1),'Tips':tip.mode}))
        dispense(merge_dicts(ethanol.liquid_class_disp,ethanol.location_disp,{'DispenseVolume': ethanol.vol_reagent,'Tips':tip.mode}))
        ethanol.set_well(s)
        dummy_n+=1
        if samples > 3:
            n=next(waste_tips_iter)
            smart_unloading(n)
        else:
            smart_unloading()
            t=Smart_loading(tips_dict,s,1)
            waste_tips=[]
            waste_tips.append(t.tips_location)
            
    reagent_dict[ethanol.name]["parameters"]["Volume_plate"]=ethanol.vol_total

    dely(30 if PCR_mode else 10)

    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Ethanol_wash 2 removal (6/9)')
    ##########################################################################################################################
    ethanol.reset_original()
    waste_column=5
    ethanol_tips_iter2=iter(waste_tips)


    for s in list_samples:
        tips_waste= next(ethanol_tips_iter2)
        load_tips(tips_waste)
        for n in range (2):
            aspirate(merge_dicts(ethanol.location_disp,ethanol.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume':round(ethanol.vol_reagent_extra,1) if n==0 else 20,'Tips':tips_waste['Tips']}))
            empty(merge_dicts(Tip_touch(Plate_1_3ml),{"Module":"POS6","TipsType":0,"Tips":tips_waste['Tips'],"Col":waste_column,"Row":1,"BottomOffsetOfZ":Z_height(Plate_1_3ml,waste_volume+ethanol.vol_reagent*2),
                    "DispenseRateOfP":20,"IfDely":True,"DelySeconds":time_to_delay(0.1),"SecondRouteRate":80}))
        smart_unloading()
        ethanol.set_well(s)
        waste_column+=1

    ethanol.reset_original()

    

    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Beads drying (7/9)')
    ##########################################################################################################################

    beads_dry=Smart_dely(180 if PCR_mode else 10)
    
    
    report(phase = 'Reagent Preparation', step = 'Elution')
    elution=mix_variables(reagent_dict,"Elution Buffer")
    tip_elution=Smart_loading(tips_dict,1,1)
    Distributing_elution=Distributing_mix(reagent_dict,elution,tip_elution,mixing=True,loops=1,empty_back=True,load_back=None)
    elution.reset_original()
    imgreplace(imgs[11])

    beads_dry.Wait()
    magnetic_down(0,0)


    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Eluting beads (8/9)')
    ##########################################################################################################################

    tips_elution=[]
    dummy_n=1 if samples==16 else 2

    for s in list_samples:
        if dummy_n>0:
            tips_elu=Smart_loading(tips_dict,s,1) 
        else:
            tips_dict=generate_tips_dict('POS4')
            tips_elu=Smart_loading(tips_dict,s,1)
        tips=tips_elu.tips_location
        mode = tips_elu.mode
        load_tips(tips)
        aspirate(merge_dicts(elution.location_disp,elution.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume': elution.vol_reagent-8,'Tips':mode}))
        dispense(merge_dicts(ethanol.location_disp,elution.liquid_class_disp,{'DispenseVolume': elution.vol_reagent-8,'Tips':mode}))
        mix(merge_dicts(ethanol.location_disp,Tip_touch(Plate_1_3ml),elution.liquid_class_mix,{"Tips":int(8),"SubMixLoopCounts":15,"MixLoopVolume":elution.vol_reagent,
                                                        "BottomOffsetOfZ":0.5,"MixOffsetOfZInLoop":Z_height(Plate_1_3ml,elution.vol_reagent),"MixOffsetOfZAfterLoop":Z_height(Plate_1_3ml,elution.vol_reagent)+1}))
        
        ethanol.set_well(s)
        elution.set_well(s)
        dummy_n-=1
        tips_elution.append(tips)
        smart_unloading(tips)

    ethanol.reset_original()

    imgreplace(imgs[12])


    dely(180 if PCR_mode else 10)
    magnetic_up(2.3,180 if PCR_mode else 10)


    ##########################################################################################################################
    report(phase = 'Clean-up', step = 'Trasnfering cleaned ligation product (9/9)')
    ##########################################################################################################################

    col=1
    tips_elution_iter=iter(tips_elution)
    elution.liquid_class_disp.update({"IfEmptyForDispense":True,"EmptyForDispenseOffsetOfZ":Z_height(PCR_plate,elution.vol_reagent-10)})
    for s in list_samples:
        tips=next(tips_elution_iter)
        mode = tips['Tips']
        load_tips(tips)
        aspirate(merge_dicts(ethanol.location_disp,elution.liquid_class_asp,{'BottomOffsetOfZ':0.5,'AspirateVolume': elution.vol_reagent-10,'Tips':mode}))
        dispense(merge_dicts(elution.liquid_class_disp,{"Module":"POS1-1","TipsType":0,"Tips":mode,"Col":col,"Row":1,'DispenseVolume': elution.vol_reagent-10,'Tips':mode}))
        ethanol.set_well(s)
        col+=1
        smart_unloading()

    imgreplace(imgs[13])

    ##########################################################################################################################
    pcr_stop_heating()
    unload_tips({'Module' : 'POS8'})
    magnetic_down(0,0)
    dialog("End of the Library Amplification\n> Purified libraries located at POS1 — store at 4°C or at -20°C.\n> Dispose of used consumables appropriately.\n> Recycle only used tips (never the box)")


##########################################################################################################################
#Variables
##########################################################################################################################

tips250ul='TipGEBAF250A'
PCR_plate='PCRBioRadHSP9601'
Plate_1_3ml='DeepwellPlateDT7350504'
Empty_box='EmptyTipGEBAF250A'
SC_micro_tube_05ml="SC_micro_tube_05ml"
SC_micro_tube_2ml="SC_micro_tube_2ml_sp100"


#PCR methods.
PCR_mode=True
PCR_start='ABL_START_25' 
PCR_25_to_4= 'ABL_25-4' if PCR_mode else PCR_start
PCR_fragmentation='ABL_MGISP100_Enzymatic_Fragmentation' if PCR_mode else PCR_start
PCR_Ligation='ABL_MGISP100_Ligation_reaction' if PCR_mode else PCR_start
PCR_Mix='ABL_MGISP100_Library_amplification_8' if PCR_mode else PCR_start


    

    
##########################################################################################################################
#Workflow
##########################################################################################################################
#workflow block begin
#Starting variables


Protocol_selection = require2({"Pop-up activated (recommended)":["Yes","No"],"Sample number":list(range(1,17)),"Select Protocol":["Adapter-Ligated Libraries","Library Amplification"]},{})
protocol=str(Protocol_selection[0][2])
log(protocol)
pop_up=True if Protocol_selection[0][0] == "Yes" else False
log(Protocol_selection[0][0])
samples = float(Protocol_selection[0][1])
log(samples)

if protocol=="Adapter-Ligated Libraries":
    library_creation(samples,PCR_mode,pop_up)
else:
    library_amplification(samples,PCR_mode,pop_up)
    


