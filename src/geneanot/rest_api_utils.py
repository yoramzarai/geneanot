# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments
"""
Utils for REST API.
"""
from functools import partial
import requests

# possible request types
request_types: dict = {
    'get': requests.get,
    'post': requests.post
}

def endpoint_base(typ: str, server: str,
                  ext: str = '',
                  params: dict | None = None,
                  headers: dict | None = None,
                  data: dict | None = None) -> dict | str:
    """Base endpoint interface."""
    assert typ in request_types, f"{typ=} not supported !!"

    if not (r := request_types[typ](f"{server}{ext}", params=params, headers=headers, data=data)).ok:
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError as err:
            print(f"Error in request {request_types[typ].__name__}: {err}")
        return {}
    try:
        return r.json()
    #except:
    except (AttributeError,requests.exceptions.JSONDecodeError):
        return r.text

# base get and post endpoints
endpoint_get_base = partial(endpoint_base, typ='get')
endpoint_post_base = partial(endpoint_base, typ='post')
