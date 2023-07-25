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

query_logger = logging.getLogger('download_log.txt')

def query_waveforms_for_agency(Agency, net, stat, channel, start, end):
    '''Query waveforms off the BGS EIDA server'''

    st = Agency.get_waveforms(network=net, station=stat, channel=channel,
                            location='00',
                            starttime=start,
                            endtime=end)
    
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
    st_downloaded = False
    n_agency = 0
    for agency in agency_list:
        client = Client(agencies[agency])
        try:
            st = query_waveforms_for_agency(client, network_code, station_code, channel,
                                            loc_code, start_time, end_time)
        except FDSNNoDataException:
            query_logger.error(f'No data available from {agency} for station {network_code}.{station_code}.{loc_code}.{channel}')
            query_logger.error(f'Event time {print(start_time)}')
            continue

        if len(st) == 3:
            # assuming we want 3 component data
            query_logger.info(f'Data found for {network_code}.{station_code}.{loc_code}.{channel} at Event time {print(start_time)}')
            # data downloaded so break loop
            return st
    

    
