import json
import requests
import os
import time
import signal

# The endpoint is static for VMs with Virtual Network.
# We should ensure our VMs have VNET.
# Discovering the endpoint for non-VNET machine does not seem trivial AND
# the example in azure docs links to an archived (since 16/11/23) git repo.
SERVICE_ENDPOINT = "http://169.254.169.254/metadata/scheduledevents"
QUERY_HEADER = {"Metadata" : "true"}
QUERY_PARAMS = {"api-version": "2020-07-01"}

SPOT_EVICTION_SIGNAL = "Preempt"

def get_scheduled_events():
    response = requests.get(SERVICE_ENDPOINT, headers=QUERY_HEADER, params=QUERY_PARAMS)
    data = response.json()
    return data["Events"]

def confirm_scheduled_event(event_id):
    payload = json.dumps({"StartRequests": [{"EventId": event_id}]})
    response = requests.post(SERVICE_ENDPOINT, headers=QUERY_HEADER, params=QUERY_PARAMS, data=payload)

def is_preempt_event(event):
    return event["EventType"] == SPOT_EVICTION_SIGNAL

def monitor_loop():
    while True:
       events = get_scheduled_events()
       for event in events:
           if is_preempt_event(event):
               return event["EventId"]

       time.sleep(2)

def kill_loop(pid):
    if os.path.exists(f"/proc/{pid}"):
        os.kill(pid, signal.SIGINT)
        while os.path.exists(f"/proc/{pid}"):
            pass

def main():
    pid = int(sys.argv[1])
    event = monitor_loop()
    kill_loop(pid)
    status = confirm_scheduled_event(event)


