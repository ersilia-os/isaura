"""Status route handlers."""


def get(event, context):
    return {
        "statusCode": 200,
        "headers": {"Content-Type": "text/plain"},
        "body": f"Hello, I am alive! You have hit........\nevent: {event}\ncontext: {str(context)}",
    }
