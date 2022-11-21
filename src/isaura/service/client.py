"""Isaura clients."""

from typing import List, Optional

import requests
from pydantic import parse_obj_as


from isaura.routes.precalc.get.app.schema import (
    ResponseBodySchema,
    QueryParams,
    QueryType,
)
from isaura.routes.schemas.common import Precalc
from isaura.utils.aws import get_dynamo_table


class IsauraAdminClient:
    def __init__(self: "IsauraAdminClient", table_name: str = "isaura_table") -> None:
        self.dynamo_table_client = get_dynamo_table(table_name)

    def insert(self: "IsauraAdminClient", precalc_list: List[Precalc]) -> None:
        """Add precals to dynamo table in batches.

        Args:
            precalc_list (List[Precalc]): List of Precalc objects
        """
        with self.dynamo_table_client.batch_writer() as batch:
            for precalc in precalc_list:
                batch.put_item(Item=precalc.dict(by_alias=True))

    def delete(self: "IsauraAdminClient", precalc_id_list: List[str]) -> None:
        """Remove precals from dynamo table in batches.

        Args:
            precalc_id_list (List[Precalc]): List of Precalc Ids to delete.
        """
        with self.dynamo_table_client.batch_writer() as batch:
            for precalc_id in precalc_id_list:
                batch.delete_item(
                    Key={"pk": precalc_id.split("#")[0], "sk": precalc_id.split("#")[1]}
                )


class IsauraClient:
    def __init__(self: "IsauraClient", api_url: str) -> None:
        self.api_url = api_url

    def get_all_precals(
        self: "IsauraClient", last_eval_key: Optional[str] = None
    ) -> ResponseBodySchema:
        """Get all precals.

        Args:
            last_eval_key (Optional[str], optional): Last eval key for pagination. Defaults to None.

        Returns:
            ResponseBodySchema: Api response.
        """
        params = QueryParams(
            query_type=QueryType.GET_ALL_PRECALC, last_eval_key=last_eval_key
        )

        resp = requests.get(url=self.api_url, params=params.dict())
        return ResponseBodySchema(
            msg=resp.json()["msg"],
            items=parse_obj_as(List[Precalc], resp.json()["items"]),
            last_eval_key=resp.json()["last_eval_key"],
            errors=resp.json()["errors"],
        )
