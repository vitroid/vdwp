from functools import wraps
import warnings
import inspect


def deprecated_alias(new_name):
    """
    名前がかわった関数を実行するデコレーター。
    """

    def decorator(func):
        # 新しい関数を取得
        frame = inspect.currentframe()
        try:
            # 呼び出し元のフレームを取得
            caller_frame = frame.f_back
            # 新しい関数を取得
            new_func = caller_frame.f_locals[new_name]
        finally:
            del frame

        @wraps(func)
        def wrapper(*args, **kwargs):
            warnings.warn(
                f"Function '{func.__name__}' is deprecated. Please use '{new_name}' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            return new_func(*args, **kwargs)

        return wrapper

    return decorator


def deprecated(new_func):
    """
    サポートが終了した関数に代わる新しい関数を推奨するデコレーター。
    """

    def decorator(old_func):
        @wraps(old_func)
        def wrapper(*args, **kwargs):
            warnings.warn(
                f"関数 '{old_func.__name__}' は非推奨です。代わりに新しい関数 '{new_func.__name__}' を使用してください。",
                DeprecationWarning,
                stacklevel=2,
            )
            return old_func(*args, **kwargs)

        return wrapper

    return decorator
