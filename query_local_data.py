import datetime
import obspy
import pandas as pd
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy import UTCDateTime
from pathlib import Path
import logging

agencies = {'BGS': 'http://eida.bgs.ac.uk', 'ORPHEUS': 'http://www.orfeus-eu.org',
            'GEOFON': 'http://geofon.gfz-potsdam.de', 'EIDA': 'http://eida-federator.ethz.ch'}

query_logger = logging.getLogger('data_download.log')

def query_waveforms_for_agency(Agency, net, stat, channel, loc_code, start, end):
    '''Query waveforms off the BGS EIDA server'''

    st = Agency.get_waveforms(network=net, station=stat, channel=channel,
                            location=loc_code,
                            starttime=start,
                            endtime=end,
                            attach_response=True)
    
    return st 

def query_eida_stat(net, stat, channel, start, end):

    bgs_eida = Client(agencies['BGS'])
    inventory = bgs_eida.get_stations(network=net, station=stat,
                                  channel=channel,
                                  level='response',
                                  starttime=UTCDateTime(start),
                                  endtime=UTCDateTime(end))
    return inventory

def multi_agency_query(agency_list, network_code, station_code, channel,
                       loc_code, start_time, end_time):
    '''
    Query waveform data by trying multiple agencies. 
    
    Parameters:
    -----------
    network_code:
        Network code of station
    station_code:
        station code
    channel:
        channel (or channels with wildcard) to query

    '''
    for agency in agency_list:
        client = Client(agencies[agency])
        try:
            st = query_waveforms_for_agency(client, network_code, station_code, channel,
                                            loc_code, start_time, end_time)
        except FDSNNoDataException:
            query_logger.warning(f'No data available from {agency} for station {network_code}.{station_code}.{loc_code}.{channel}. Event time {start_time}')
            continue

        if len(st) == 3:
            # assuming we want 3 component data
            query_logger.info(f'Data found for {network_code}.{station_code}.{loc_code}.{channel} at Event time {start_time}')
            # data downloaded so break loop
            return st
        elif len(st) % 3 ==0:
            #In some cases we get duplicates
            query_logger.warning(f'Duplicate (possibly) channels returned for for {network_code}.{station_code} for event at {start_time}. Only returning first set.')
            chs = [channel.strip('?') + cmp for cmp in ['Z','N','E']]
            trZ = st.select(channel=chs[0])[0]
            trN = st.select(channel=chs[1])[0]
            trE = st.select(channel=chs[2])[0]
            return obspy.Stream([trZ, trN, trE])
        else:
            query_logger.warning(f'We do not get 3-component data from {agency} for {network_code}.{station_code} for event at {start_time}')

    return 'No Data'
    
