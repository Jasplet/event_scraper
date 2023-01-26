
import argparse
from datetime import timedelta
import obspy
from obspy.clients.fdsn.client import Client
from obspy.clients.fdsn.header import (FDSNNoDataException)
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth, kilometer2degrees
from obspy.io.sac import SACTrace
from obspy.taup import TauPyModel
import warnings
from os.path import isfile

skipped = 0
# Query ISC catalog (via obspy) for a set of events a set distance from
# a station then using this Catalog object, query waveforms where we would
# expect to see SKS.

def make_event_query(station, t1, t2, phase='SKS'):
    '''
    Function to query events which should give a clear 
    teleseismic phase. Currently only implemented for SKS
    '''
    client = Client('IRIS')
    if phase == 'SKS':
        event_catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=5.5,
        maxmagnitude=7, catalog="ISC", latitude=station['latitude'], longitude=station['longitude'],
        minradius=95, maxradius=120, magnitudetype="MW")

    return event_catalog

def query_waveforms(Event, station, loc_code="00"):
    ''' 
    Use Obspy Client.get_waveforms() to query waveform data.
    Quries waveforms up to 30 minutes after origin time from a single 
    earthquake (Obspy Event object)
    '''

    client = Client("IRIS")
    # Need to make a list and do individual queries for a single station
    # as we cannot easily split up the Stream object from get_waveforms_bulk 
    # so insteasd we will make each Event's query seperately and compile
    # a list of (in theory 3 component) Stream objects
    warnings.filterwarnings('error')
    evdp = Event.origins[0].depth / 1000 # converts from [m] -> [km]
    dist_m, _, _ = gps2dist_azimuth(Event.origins[0].latitude, Event.origins[0].longitude, station['latitude'], station['longitude'])
    dist_deg = kilometer2degrees(dist_m/1000)
    tt_pred = calc_tt(evdp, dist_deg)
    # Request records that start 3 minutes before and and 3 minutes after the predicted SKS arrival
    # SKS arrival is calculated in seconds after 
    t1 = Event.origins[0].time + timedelta(seconds=tt_pred) - timedelta(minutes=3)
    t2 = Event.origins[0].time + timedelta(seconds=tt_pred) + timedelta(minutes=3)
    st_raw = client.get_waveforms("IU", station['name'], loc_code, "BH?",t1, t2,
                                    minimumlength=360, attach_response=True)
    if len(st_raw) > 3:
        raise ValueError('Too many channels')
    channels = [tr.stats.channel for tr in st_raw]
    if ('BH1' in channels) and ('BH2' in channels):
        # get invetory fix orientations
        inv = client.get_stations(network='IU', sta=station['name'], loc=loc_code, channel="BH?",
                                starttime=t1, endtime=t2, level='response') 
        st_raw.rotate('->ZNE', inventory=inv)
    #Copy to avoid accidental double correcting
    try:
        st = st_raw.remove_response().copy() #remove instrument reponse now (might as well)
    except ValueError:
        print('no reponse, so skip event')
        return   
    # Intial broad fliter, we hont want any frequency information outside this range 
    # (and will likely want to get rid of some of the high f noise (0.1/0.3-0.5Hz))
    st.filter('bandpass', freqmin=0.01, freqmax=0.5, zerophase=True)

    return st 

def get_waveforms_for_catalog(Catalog, station, loc_code="00", write_out=True, path=None):
    '''
    Queries waveform data for a Catalog of Events at a single station.
    Waveform data can be written out (as SAC files) or
    returned as a single Stream
    '''
    print(f'{len(Catalog)} events to query')
    # turn UserWarnings (which we get if there is no response into errors
    # which we can try, except. 
    if write_out:
        #To write out the Stream objects as sac files with updated headers we need to
        #do a bit of a hack
        skipped = 0
        n = 0
        for Event in Catalog:
            origin = Event.origins[0].time
            filename = f'{station["name"]}_{origin.year}{origin.julday:03d}_{origin.hour:02d}{origin.minute:02d}{origin.second:02d}'
            if isfile(f'{path}/{filename}.BHN') & isfile(f'{path}/{filename}.BHE') & isfile(f'{path}/{filename}.BHZ'):
                print('{filename} already downloaded')
                n += 1
                continue
            else:
                try:
                    st = query_waveforms(Event, station, loc_code)
                
                    for tr in st:
                        channel = tr.stats.channel
                        sactr = make_sactrace(Event, station, tr)
                        
                        sactr.write(f'{path}/{filename}.{channel}', byteorder='big')
                    n +=1 
                    if n % 10 == 0:
                        print(n)
                except UserWarning:
                    print('f{UserWarning}, skip event')
                    skipped +=1
                    continue
                except FDSNNoDataException:
                    print(f'{FDSNNoDataException}, skip event')
                    skipped +=1
                    continue

        print(f'{skipped}/{len(Catalog)} event skipped for having no responses')
        return 
    else:
        st_out = obspy.Stream()
        for Event in Catalog:
            st = query_waveforms(Event, station, loc_code)
            st_out += st
        return st_out
            

def make_sactrace(Event, station, Trace):
    '''
    use the Obspy SACTrace object to update sac headers more accurately
    '''
    sactr = SACTrace.from_obspy_trace(Trace)
    # Set origin time
    otime = Event.origins[0].time
    evla = Event.origins[0].latitude
    evlo = Event.origins[0].longitude
    evdp = Event.origins[0].depth / 1000 # converts from [m] -> [km]
    dist_m, az, baz = gps2dist_azimuth(evla, evlo, station['latitude'], station['longitude'])
    dist_deg = kilometer2degrees(dist_m/1000)
    sks_rel_tt = calc_tt(evdp, dist_deg)
    sactr.o = otime
    sactr.iztype = 'io'
    sactr.evla = evla
    sactr.evlo = evlo
    sactr.evdp = evdp
    sactr.stla = station['latitude']
    sactr.stlo = station['longitude']
    sactr.kstm = station['name']
    sactr.gcarc = dist_deg
    sactr.az = az
    sactr.baz = baz
    # Add some initial windows based off 
    # 1-D predicted SKS traveltime
    sactr.user0 = sks_rel_tt - 2.5
    sactr.user1 = sks_rel_tt + 2.5
    sactr.user2 = sks_rel_tt + 12.5
    sactr.user3 = sks_rel_tt + 17.5

    return sactr

def calc_tt(evdp, dist_deg, phase='SKS'):

    model = TauPyModel('iasp91')
    arrivals = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=dist_deg, phase_list=[phase])
    if len(arrivals) > 1:
        raise ValueError('2 SKS phases!!!')
    return arrivals[0].time 

if __name__ == '__main__':

    # parser = argparse.ArgumentParser(
    #     description="Queries ISC bulletin (via IRIS) to get events that should give SKS phases at a station")
    # parser.add_argument('-s','-station',action='store',type=str, help='Station Code')
    # parser.add_argument('-la', '-latitude',action='store',type=float,help='station latitude')
    # parser.add_argument('-lo', '-longitude',action='store',type=float,help='station longitude')
    station = {'name':'FURI', 'latitude': 8.895, 'longitude': 38.86} 
    start = UTCDateTime("2001-01-01T00:00:00")
    end = UTCDateTime("2022-01-01T00:00:00")

    path = '/Users/ja17375/Projects/DeepMelt/Ethiopia/FURI_data/'
    # events = make_event_query(station, start, end)
    # events.write(f"{path}/{station}_data.xml", format="QUAKEML")
    events = obspy.read_events(f'{path}/FURI_data.xml', format='QUAKEML')
    get_waveforms_for_catalog(events, station, write_out=True,
                path='/Users/ja17375/Projects/DeepMelt/Ethiopia/FURI_data/data')
   
