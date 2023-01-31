
import argparse
from datetime import timedelta
import obspy
from obspy.clients.fdsn.client import Client
from obspy import UTCDateTime
from obspy.geodetics import locations2degrees


# Query ISC catalog (via obspy) for a set of events a set distance from
# a station then using this Catalog object, query waveforms where we would
# expect to see SKS.

def make_event_query(stla, stlo, t1, t2, phase='SKS'):

    client = Client('IRIS')
    if phase == 'SKS':
        event_catalog = client.get_events(starttime=t1, endtime=t2, minmagnitude=5.5,
        maxmagnitude=7, catalog="ISC", latitude=stla, longitude=stlo,
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
    t1 = event.origins[0].time
    t2 = t1 + timedelta(minutes=30)
    st = client.get_waveforms("IU", station, loc_code, "BH?",t1, t2,
                                minimumlength=1800, attach_response=True)
    return st 

def get_waveforms_for_catalog(Catalog, station, loc_code="00", write_out=True, path=None):
    '''
    Queries waveform data for a Catalog of Events at a single station.
    Waveform data can be written out (as SAC files) or
    returned as a single Stream
    '''

    for Event in Catalog:
        st = query_waveforms(Event, station, loc_code)
        if write_out:
            #To write out the Stream objects as sac files with updated headers we need to
            #do a bit of a hack


def write_out_st():
    '''
    Writes out the list of Streams to individual SAC files
    '''

def update_sac_headers():
    '''
    Update sac headers for each Stream object
    '''


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Queries ISC bulletin (via IRIS) to get events that should give SKS phases at a station")
    parser.add_argument('-s','-station',action='store',type=str, help='Station Code')
    parser.add_argument('-la', '-latitude',action='store',type=float,help='station latitude')
    parser.add_argument('-lo', '-longitude',action='store',type=float,help='station longitude')
    station = 'FURI' # add argparse later
    stla = '8.895'
    stlo = '38.68'
    start = UTCDateTime("2001-01-01T00:00:00")
    end = UTCDateTime("2022-01-01T00:00:00")

    path = '/Users/ja17375/Projects/DeepMelt/Ethiopia/FURI_data/'
    events = make_event_query(stla, stlo, start, end)
    events.write(f"{path}/{station}_data.xml", format="QUAKEML")

